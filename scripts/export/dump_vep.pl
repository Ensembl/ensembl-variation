#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.



=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk.org>.

=cut


use Getopt::Long;
use DBI;

my $config = {};

my %special_options = (
  'homo_sapiens'      => ' --sift b --polyphen b --regulatory'.
                         ' --freq_file /nfs/ensembl/wm2/VEP/cache/ALL_1KG_ESP_freqs_with_alleles.txt,AFR,AMR,ASN,EUR,AA,EA',
  'mus_musculus'      => ' --regulatory --sift b',
  'bos_taurus'        => ' --sift b',
  'canis_familiaris'  => ' --sift b',
  'danio_rerio'       => ' --sift b',
  'equus_caballus'    => ' --sift b',
  'gallus_gallus'     => ' --sift b',
  'rattus_norvegicus' => ' --sift b',
  'sus_scrofa'        => ' --sift b',
  'triticum_aestivum' => ' --sift b',
);

GetOptions(
    $config,
  'hosts|h=s',
  'user|u=s',
  'port|P=i',
  'password|p=s',
  'dir=s',
  'command=s',
  'mem=i',
  'pattern=s',
  'exclude_pattern=s',
  'queue=s',
  'version=i',
  'overwrite',
  'refseq',
) or die "ERROR: Could not parse command line options\n";

# set defaults
$config->{hosts}    ||= 'ens-staging,ens-staging2';
$config->{user}     ||= 'ensro';
$config->{port}     ||= 3306;
$config->{dir}      ||= $ENV{'HOME'}.'/.vep/';
$config->{command}  ||= 'perl variant_effect_predictor.pl --build all';
$config->{mem}      ||= 12000;
$config->{queue}    ||= 'normal';

# check dir exists
die "ERROR: Dump directory ".$config->{dir}." does not exist\n" unless -e $config->{dir};

# check version defined
die "ERROR: No Ensembl DB version defined - use --version [version]\n" unless defined($config->{version});

# check command supplied looks sensible
die "ERROR: Supplied command doesn't look right, it should look something like:\n\nperl -I /include/perl/libs/ /path/to/variant_effect_predictor.pl --build all\n\n"
  unless $config->{command} =~ /variant_effect_predictor.+\-build (all|\w+)/;

# parse hosts
$config->{hosts} = [split /\,/, $config->{hosts}];

foreach my $host(@{$config->{hosts}}) {
  debug("Getting species list for host $host");
  my $species_list = get_species_list($config, $host);
  debug("Found ".(scalar @$species_list)." valid databases\n");
  
  foreach my $species_hash(@$species_list) {
    my ($species, $assembly) = ($species_hash->{species}, $species_hash->{assembly});
    
    debug("Running VEP dump for $species");
    dump_vep($config, $host, $species, $assembly);
    
    debug("Compressing dump file");
    tar($config, $species, $assembly);
  }
}

sub get_species_list {
  my $config = shift;
  my $host   = shift;

  my $connection_string = sprintf(
      "DBI:mysql(RaiseError=>1):host=%s;port=%s",
      $host,
      $config->{port}
    );
  
  # connect to DB
  my $dbc = DBI->connect(
      $connection_string, $config->{user}, $config->{password}
  );
  
  my $version = $config->{version};

  my $sth = $dbc->prepare(qq{
    SHOW DATABASES LIKE '%\_core\_$version%'
  });
  $sth->execute();
  
  my $db;
  $sth->bind_columns(\$db);
  
  my @dbs;
  push @dbs, $db while $sth->fetch;
  $sth->finish;
  
  # refseq?
  if(defined($config->{refseq})) {
    $sth = $dbc->prepare(qq{
      SHOW DATABASES LIKE '%\_otherfeatures\_$version%'
    });
    $sth->execute();
    $sth->bind_columns(\$db);
    
    push @dbs, $db while $sth->fetch;
    $sth->finish;
  }

  # remove master and coreexpression
  @dbs = grep {$_ !~ /master|express/} @dbs;

  # filter on pattern if given
  my $pattern = $config->{pattern};
  my $exclude = $config->{exclude_pattern};
  @dbs = grep {$_ =~ /$pattern/} @dbs if defined($pattern);
  @dbs = grep {$_ !~ /$exclude/} @dbs if defined($exclude);

  my @species;

  foreach my $current_db_name (@dbs) {
    
    # special case otherfeatures
    # check it has refseq transcripts
    if($current_db_name =~ /otherfeatures/) {
    
      # get assembly name
      $sth = $dbc->prepare("select meta_value from ".$current_db_name.".meta where meta_key='assembly.default';");
      $sth->execute();
      my $assembly = $sth->fetchall_arrayref()->[0]->[0];
      die("ERROR: Could not get assembly name from meta table for $current_db_name\n") unless $assembly;
      
      $sth = $dbc->prepare(qq{
        SELECT COUNT(*)
        FROM $current_db_name\.transcript
        WHERE stable_id LIKE 'NM%'
        OR source = 'refseq'
      });
      $sth->execute;
      
      my $count;
      $sth->bind_columns(\$count);
      $sth->fetch;
      
      if($count) {
        $current_db_name =~ s/^([a-z]+\_[a-z,1-9]+)(\_[a-z]+)?(.+)/$1$2/;
        $current_db_name =~ s/\_otherfeatures$/\_refseq/;
        push @species, { species => $current_db_name, assembly => $assembly};
      }
    }
    
    else {
      # get assembly and species names
      $sth = $dbc->prepare("select species_id, meta_key, meta_value from ".$current_db_name.".meta where meta_key in ('assembly.default', 'species.production_name');");
      $sth->execute();
      
      my ($species_id, $key, $value, $by_species);
      $sth->bind_columns(\$species_id, \$key, \$value);
      
      $by_species->{$species_id}->{$key} = $value while $sth->fetch();
      $sth->finish();
      
      my $count = 0;
      
      foreach my $hash(values %$by_species) {
        next unless $hash->{'assembly.default'} && $hash->{'species.production_name'};
        push @species, { species => $hash->{'species.production_name'}, assembly => $hash->{'assembly.default'}};
        $count++;
      }
      
      die("ERROR: Problem getting species and assembly names from $current_db_name; check meta table\n") unless $count;
    }
  }

  return \@species;
}

sub dump_vep {
  my $config = shift;
  my $host   = shift;
  my $sp     = shift;
  my $ass    = shift;
  
  # check if dir exists
  if(!defined($config->{overwrite}) && -e $config->{dir}.'/'.$sp.'/'.$config->{version}.'_'.$ass) {
    debug("Existing dump directory found for $sp $ass, skipping (use --overwrite to overwrite)\n");
    return;
  }
  
  my $species = $sp;
  my $refseq = $species =~ s/\_refseq$//;
  
  my $command = join " ", (
    sprintf(
      'bsub -K -J %s_%s -M %i -R"select[mem>%i] rusage[mem=%i]" -q %s -o %s -e %s',
      $species,
      $ass,
      $config->{mem},
      $config->{mem},
      $config->{mem},
      $config->{queue},
      $config->{dir}.'/'.$sp.'_'.$ass.'_vep_dump.farmout',
      $config->{dir}.'/'.$sp.'_'.$ass.'_vep_dump.farmerr',
    ),
    $config->{command},
    defined $special_options{$species} ? $special_options{$species} : "",
    "--no_adaptor",
    "--species ".$species,
    "--host ".$host,
    "--user ".$config->{user},
    "--port ".$config->{port},
    "--dir ".$config->{dir},
    defined $config->{password} ? "--password ".$config->{password} : "",
    $refseq ? '--refseq' : '',
  );
  
  debug("Use \"tail -f ".$config->{dir}.'/'.$species.'_vep_dump.farmout'."\" to check progress");
  system($command);
  #print "$command\n";
}

sub tar {
  my $config   = shift;
  my $species  = shift;
  my $assembly = shift;
  
  my $tar_file = $config->{dir}."/".$species."_vep_".$config->{version}."_".$assembly.".tar.gz";
  
  # check if tar exists
  if(!defined($config->{overwrite}) && -e $tar_file) {
    debug("Existing dump file found for $species, skipping (use --overwrite to overwrite)\n");
    return;
  }
  
  # check dir exists
  my $root_dir = $config->{dir};
  my $sub_dir  = $species."/".$config->{version}."_".$assembly;
  
  die("ERROR: VEP dump directory $root_dir/$sub_dir not found") unless -e $root_dir.'/'.$sub_dir;
  
  my $command = "tar -cz -C $root_dir -f $tar_file $sub_dir";
  system($command);# or die "ERROR: Failed to create tar file $tar_file\n";
  #print "$command\n";
}

# gets time
sub get_time() {
  my @time = localtime(time());

  # increment the month (Jan = 0)
  $time[4]++;

  # add leading zeroes as required
  for my $i(0..4) {
    $time[$i] = "0".$time[$i] if $time[$i] < 10;
  }

  # put the components together in a string
  my $time =
    ($time[5] + 1900)."-".
    $time[4]."-".
    $time[3]." ".
    $time[2].":".
    $time[1].":".
    $time[0];

  return $time;
}

# prints debug output with time
sub debug {
  my $text = (@_ ? (join "", @_) : "No message");
  my $time = get_time;
  
  print $time." - ".$text.($text =~ /\n$/ ? "" : "\n");
}
