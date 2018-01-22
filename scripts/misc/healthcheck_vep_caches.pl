#!/usr/bin/env perl

=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

use strict;
use Getopt::Long;
use Storable qw(fd_retrieve);
use Bio::EnsEMBL::Registry;
use Data::Dumper;
use IPC::Open3;

$Data::Dumper::Maxdepth = 1;
$Data::Dumper::Indent = 1;

BEGIN {
  $| = 1;
  use Test::More;
}

my $config = {};

GetOptions(
  $config,
  'help',                  # displays help message
  
  'dir|d=s',                 
  'version|v=i',
  'cache_version|c=i',
  'species|s=s',
  
  'host|h=s',                  # DB options
  'user|u=s',
  'password=s',
  'port|p=i',
  
  'no_fasta|nf',   # don't check FASTA files
  
  'random|r=f',
  'max_vars|m=i',
) or die "ERROR: Failed to parse command-line flags\n";

do { usage(); exit(0); } if $config->{help};

$config->{reg} = 'Bio::EnsEMBL::Registry';
$config->{reg}->no_version_check(1);

$config->{dir}           ||= $ENV{HC_VEP_DIR}           || $ENV{HOME}.'/.vep/';
$config->{version}       ||= $ENV{HC_VEP_VERSION}       || $config->{reg}->software_version;
$config->{cache_version} ||= $ENV{HC_VEP_CACHE_VERSION} || $config->{version};
$config->{species}       ||= $ENV{HC_VEP_SPECIES}       || 'all';
$config->{host}          ||= $ENV{HC_VEP_HOST}          || 'ens-staging1,ens-staging2';
$config->{user}          ||= $ENV{HC_VEP_USER}          || 'ensro';
$config->{port}          ||= $ENV{HC_VEP_PORT}          || 3306;
$config->{password}      ||= $ENV{HC_VEP_PASS}          || undef;
$config->{no_fasta}      ||= $ENV{HC_VEP_NO_FASTA}      || undef;
$config->{max_vars}      ||= $ENV{HC_VEP_MAX_VARS}      || 100;

if(!defined($config->{random})) {
  if(defined($ENV{HC_VEP_RANDOM})) {
    $config->{random} = $ENV{HC_VEP_RANDOM};
  }
  else {
    $config->{random} = 0;
  }
}

$config->{cache_region_size} = 1e6;


if(defined($config->{random})) {
  die("ERROR: --random must be between 0 and 1\n") unless $config->{random} >= 0 && $config->{random} <= 1;
}

my %match_species = ();
if($config->{species} ne 'all') {
  $match_species{$_} = 1 for split(/\,/, $config->{species});
}

foreach my $host(split /\,/, $config->{host}) {
  
  my $species_list = get_species_list($config, $host);
  
  if($config->{species} eq 'all') {
    $match_species{$_} = 1 for map {$_->{species}} @$species_list;
  }
  
  $config->{reg}->load_registry_from_db(
    -host       => $host,
    -user       => $config->{user},
    -pass       => $config->{password},
    -port       => $config->{port},
    -db_version => $config->{version},
  );
  
  foreach my $species_hash(@$species_list) {
    
    my ($species, $assembly) = ($species_hash->{species}, $species_hash->{assembly});
    
    next unless $match_species{$species};
    
    print STDERR "\nChecking $species\n";
    
    $config->{tests}->{$species} = {};
    $config->{current} = $config->{tests}->{$species};
    
    my $sp = $species;
    my $refseq = $sp =~ s/\_refseq//;
    
    my $has_var = has_variation($config, $sp);
    my $has_reg = has_regulation($config, $sp);
    
    my $dir = $config->{dir}.'/'.$species.'/'.$config->{cache_version}.($config->{version} >= 76 ? '_'.$assembly : '');
    
    ok(-d $config->{dir}.'/'.$species, "\[$species\] species dir exists");
    ok(-d $dir, "\[$species\] version dir exists");
    
    ok(-e $dir.'/info.txt', "\[$species\] info.txt exists");
    ok(open(IN, $dir.'/info.txt'), "\[$species\] info.txt readable");
    
    my ($v, $r, $tabix);
    
    while(<IN>) {
      ok(/all/, "--build all used") if /^build/;
      $v = $_ if /^variation_cols/;
      $r = $_ if /^regulatory/;
      $tabix = $_ if /^var_type/;
      $config->{current}->{expect_sift} = 1 if /sift_version/i;
      $config->{current}->{expect_polyphen} = 1 if /polyphen_version/i;
    }
    close IN;
    
    ok(defined($v), "\[$species\] variation_cols defined in info.txt") if $has_var;
    ok(defined($r), "\[$species\] regulation defined in info.txt") if $has_reg;
    
    # store var cols etc
    if($has_var && $v) {
      $config->{current}->{variation_cols} = [split(",", (split("\t", $v))[-1])];
      $config->{current}->{vdb} = $config->{reg}->get_adaptor($sp, 'variation', 'variation')->db->dbc;
    }
        
    my @types = ('');
    push @types, '_var' if $has_var && !$tabix;
    push @types, '_reg' if $has_reg;
    
    my $sa = $config->{reg}->get_adaptor($sp, $refseq ? 'otherfeatures' : 'core', 'slice');
    my @slices = @{$sa->fetch_all('toplevel')};
    
    my @missing_dirs = ();
    my @missing_files = ();
    my @tabix_broken = ();
    
    my $i = 0;
    
    # get transcript counts
    my $sth = $sa->{dbc}->prepare(qq{
      SELECT seq_region_id, COUNT(*)
      FROM transcript
      GROUP BY seq_region_id;
    });
    $sth->execute();
    
    my ($sr_id, $count, %counts);
    $sth->bind_columns(\$sr_id, \$count);
    $counts{$sr_id} += $count while $sth->fetch;
    $sth->finish;

    # get varfeat counts
    if($has_var && $v) {
      $sth = $config->{current}->{vdb}->prepare(qq{
        SELECT seq_region_id, COUNT(*)
        FROM variation_feature
        GROUP BY seq_region_id;
      });
      $sth->execute();

      $sth->bind_columns(\$sr_id, \$count);
      $counts{$sr_id} += $count while $sth->fetch;
      $sth->finish;
    }
    
    foreach my $slice(@slices) {
      # printf("\r - checking chromosome dirs %i / %i", ++$i, scalar @slices);
      
      my $chr = $slice->seq_region_name;
      
      next unless $counts{$slice->get_seq_region_id};
      
      if(!-d $dir.'/'.$chr) {
        push @missing_dirs, $chr;
        next;
      }
      
      # initial region
      my $start = 1 + ($config->{cache_region_size} * int($slice->start / $config->{cache_region_size}));
      my $end   = ($start - 1) + $config->{cache_region_size};
      
      my $region_count = int($config->{random} * $slice->end / $config->{cache_region_size}) + 1;
      my $checked_count = 0;
      
      while($start < $slice->end) {
        
        foreach my $type(@types) {
          my $file = $dir.'/'.$chr.'/'.$start.'-'.$end.$type.'.gz';
          push @missing_files, $file unless -e $file;
        }
        
        # check this region?
        if($config->{random} && ($checked_count < $region_count) && (rand() < $config->{random})) {

          check_region(\@types, $dir, $chr, $start, $end);

          $checked_count++;
        }
        
        # increment by cache_region_size to get next region
        $start += $config->{cache_region_size};
        $end += $config->{cache_region_size};
      }
     
      # check tabix
      if($tabix) {
        push @tabix_broken, $chr unless -e $dir.'/'.$chr.'/all_vars.gz';
        push @tabix_broken, $chr.' (index)' unless -e $dir.'/'.$chr.'/all_vars.gz.tbi';
        
        if(`which tabix 2>&1` =~ /tabix$/) {
          my $filepath = $dir.'/'.$chr.'/all_vars.gz';
          my $result = `tabix $filepath 1:1-1 2>&1`;
          push @tabix_broken, $chr.' ('.$result.')' if $result =~ /fail/;
          
          if($config->{random}) {
            my $vs = 0;
            
            my $pid = open3(\*OUT, \*IN, \*ERR, "gzip -qdc $filepath");
            while(<IN>) {
              next unless rand() < $config->{random};
              chomp;
              check_variation(parse_variation($config, $_));
              if(++$vs > $config->{max_vars}) {
                kill -9, $pid;
                last;
              }
            }
            close IN;
            close OUT;
            close ERR;
          }
        }
      }
    }
    
    # only want to display a few if loads
    @missing_dirs = (@missing_dirs[0..4], "...and ".(scalar @missing_dirs - 5)." more...") if scalar @missing_dirs > 5;
    @missing_files = (@missing_files[0..4], "...and ".(scalar @missing_files - 5)." more...") if scalar @missing_files > 5;
    @tabix_broken = (@tabix_broken[0..4], "...and ".(scalar @tabix_broken - 5)." more...") if scalar @tabix_broken > 5;
    
    ok(scalar @missing_dirs == 0, "\[$species\] toplevel slice directories present")
      or diag("Missing directories for the following chromosome names:\n".join(", ", @missing_dirs)."\n");
    ok(scalar @missing_files == 0, "\[$species\] cache files present")
      or diag("Missing the following cache files:\n".join(", ", @missing_files)."\n");
    
    
    if($tabix) {
      ok(scalar @tabix_broken == 0, "\[$species\] tabix file checks")
        or diag("Tabix-based variation cache broken for the following chromosome names:\n".join(", ", @tabix_broken)."\n");
    }
    
    # sift/polyphen
    if($config->{current}->{found_protein_coding}) {
      foreach my $tool(grep {$config->{current}->{'expect_'.$_}} qw(sift polyphen)) {
        ok($config->{current}->{'found_'.$tool}, "found $tool data");
      }
    }
    
    
    # check FASTA file
    unless($config->{no_fasta}) {
      opendir CACHE, $dir;
      my @read = readdir CACHE;
      my ($fa)  = grep {/\.fa$/} @read;
      my ($idx) = grep {/\.index$/} @read;
      closedir CACHE;
    
      ok(defined($fa), "\[$species\] FASTA file present");
      ok(defined($idx), "\[$species\] FASTA index file present");
      ok($fa =~ /^$species/i, "\[$species\] FASTA file matches species name") if($fa);
      ok($fa =~ /\.$config->{version}\./, "\[$species\] FASTA file matches version") if($idx);
    
      # only check indexing ok if built already
      if(defined($idx)) {
        eval q{
          use Bio::DB::Fasta;
        };
      
        if(!$@) {
          eval q{
            Bio::DB::Fasta->new($dir.'/'.$fa);
          };
        
          ok(!$@, "\[$species\] FASTA index loads in Bio::DB::Fasta") or diag($@);
        }
      }
    }
  }
}

done_testing();

sub get_species_list {
  my $config = shift;
  my $host   = shift;

  my $connection_string = sprintf(
    "DBI:mysql(RaiseError=>1):host=%s;port=%s",
    $host,
    $config->{port}
  );
  
  # connect to DB
  $config->{dbc} = DBI->connect(
    $connection_string, $config->{user}, $config->{password}
  );
  
  my $version = $config->{cache_version};

  my $sth = $config->{dbc}->prepare(qq{
    SHOW DATABASES LIKE '%\_core\_$version%'
  });
  $sth->execute();
  
  my $db;
  $sth->bind_columns(\$db);
  
  my @dbs;
  push @dbs, $db while $sth->fetch;
  $sth->finish;
  
  # refseq?
  $sth = $config->{dbc}->prepare(qq{
    SHOW DATABASES LIKE '%\_otherfeatures\_$version%'
  });
  $sth->execute();
  $sth->bind_columns(\$db);
  
  push @dbs, $db while $sth->fetch;
  $sth->finish;

  # remove master and coreexpression
  @dbs = grep {$_ !~ /master|express/} @dbs;

  # filter on pattern if given
  my $pattern = $config->{pattern};
  @dbs = grep {$_ =~ /$pattern/} @dbs if defined($pattern);

  my @species;
  my $dbc = $config->{dbc};

  foreach my $current_db_name (@dbs) {
    
    # special case otherfeatures
    # check it has refseq transcripts
    if($current_db_name =~ /otherfeatures/) {
      
      $sth = $dbc->prepare("SELECT version FROM ".$current_db_name.".coord_system ORDER BY rank LIMIT 1;");
      $sth->execute();
      my $assembly = $sth->fetchall_arrayref()->[0]->[0];
      die("ERROR: Could not get assembly name from meta table for $current_db_name\n") unless $assembly;
      
      $sth = $config->{dbc}->prepare(qq{
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
      $sth = $dbc->prepare("select species_id, meta_value from ".$current_db_name.".meta where meta_key = 'species.production_name';");
      $sth->execute();
      
      my ($species_id, $value, $species_ids);
      $sth->bind_columns(\$species_id, \$value);
      
      $species_ids->{$species_id} = $value while $sth->fetch();
      $sth->finish();
      
      my $count = 0;
      
      foreach $species_id(keys %$species_ids) {
        $sth = $dbc->prepare("SELECT version FROM ".$current_db_name.".coord_system WHERE species_id = ".$species_id." ORDER BY rank LIMIT 1;");
        $sth->execute();
        my $assembly;
        $sth->bind_columns(\$assembly);
        $sth->execute();
        $sth->fetch();
        
        next unless $assembly;
        
        # copy server details
        my %species_hash;
        
        $species_hash{species} = $species_ids->{$species_id};
        $species_hash{assembly} = $assembly;
        $sth->finish();

        push @species, \%species_hash;
        $count++;
      }
      
      die("ERROR: Problem getting species and assembly names from $current_db_name; check meta table\n") unless $count;
    }
  }

  return \@species;
}

sub has_variation {
  my ($config, $species) = @_;
  
  my $version = $config->{cache_version};
  
  my $sth = $config->{dbc}->prepare(qq{
    SHOW DATABASES LIKE '$species\_variation\_$version%'
  });
  $sth->execute();
  
  my $db;
  $sth->bind_columns(\$db);
  $sth->fetch;
  $sth->finish;
  
  return defined($db);
}

sub has_regulation {
  my ($config, $species) = @_;
  
  my $version = $config->{cache_version};
  
  my $sth = $config->{dbc}->prepare(qq{
    SHOW DATABASES LIKE '$species\_funcgen\_$version%'
  });
  $sth->execute();
  
  my $db;
  $sth->bind_columns(\$db);
  $sth->fetch;
  $sth->finish;
  
  return 0 unless defined($db);
  
  my $v = $config->{dbc}->selectall_arrayref(qq{
    SELECT version FROM $db.regulatory_build
  })->[0]->[0];
  
  return defined($v) && $v > 0;
}

sub check_region {
  my ($types, $dir, $chr, $start, $end) = @_;

  foreach my $type(@$types) {
    my $short_name = $chr.'/'.$start.'-'.$end.$type.'.gz';
    my $file = $dir.'/'.$short_name;
    
    if($type eq '_var') {
      my $vs = 0;
      
      my $pid = open3(\*OUT, \*IN, \*ERR, "gzip -qdc $file");
      while(<IN>) {
        next unless rand() < $config->{random};
        chomp;
        check_variation(parse_variation($config, $_));
        if(++$vs > $config->{max_vars}) {
          kill -9, $pid;
          last;
        }
      }
      close IN;
      close OUT;
      close ERR;
    }
    
    else {
      open my $fh, "gzip -dc $file |";
      my $data;
      ok($data = fd_retrieve($fh), "Open cache file $short_name") or diag($!);
      close $fh;

      # should be a hash keyed on chr name
      if(ok($data && ref($data) eq 'HASH', "Cache file is hashref")) {

        ok($data->{$chr}, "Cache file has chr key") or diag("expected $chr, got ".join(",", keys %$data));

        my $chr_data = $data->{$chr};

        # reg cache is a hash with a key for each feature type
        if($type eq '_reg') {
          ok(ref($chr_data) eq 'HASH', "Reg feat cache is hashref");

          foreach my $key(keys %$chr_data) {
            ok($key eq 'RegulatoryFeature' || $key eq 'MotifFeature', "key is RegulatoryFeature or MotifFeature") or diag("got $key");

            my $val = $chr_data->{$key};
            if(ok(ref($val) eq 'ARRAY', "value is arrayref")) {
              check_rf($_, $key) for @$val;
            }
          }
        }

        # transcript cache is an array
        else {
          if(ok(ref($chr_data) eq 'ARRAY', "Transcript cache is arrayref")) {
            check_tr($_) for @$chr_data;
          }
        }
      }
    }
  }
}

sub check_basics {
  my $ref = shift;
  
  my $n = $ref->{dbID} || $ref->{variation_name};
  
  ok($ref->{start}, "Object has start") or diag(Dumper $ref);
  ok($ref->{end}, "Object has end") or diag(Dumper $ref);
  ok($ref->{dbID} || $ref->{variation_name}, "Object has ID") or diag(Dumper $ref);
  ok(defined($ref->{strand}) || defined($ref->{seq_region_strand}), "Object has strand") or diag(Dumper $ref);
}

sub check_variation {
  my $v = shift;
  
  check_basics($v);
  
  my $n = $v->{variation_name};
  
  # check frequency
  if($v->{minor_allele}) {
    ok($v->{minor_allele} =~ /^[ACGTN-]+$/, "minor allele is valid allele") or diag("minor allele ".$v->{minor_allele});
    ok($v->{minor_allele_freq} =~ /^[0-9\.]+$/ && $v->{minor_allele_freq} >= 0 && $v->{minor_allele_freq} <= 0.5, "$n 0 <= frequency <= 0.5") or diag("$n freq ".$v->{minor_allele_freq});
  }
  
  # check population freqs
  foreach my $p(qw(AFR AMR ASN EAS EUR SAS AA EA)) {
    if($v->{$p}) {
      foreach my $set(split(',', $v->{$p})) {
        my ($a, $f) = split(':', $set);
        
        if(!defined($f)) {
          $f = $a;
          $a = undef;
        }
        
        if($a) {
          ok($a =~ /^[ACGTN-]+$/, "$p allele is valid allele") or diag("$n $p minor allele $a");
        }
        
        ok($f =~ /^[0-9\.]+$/ && $f >= 0 && $f <= 1, "0 <= frequency <= 1") or diag("$n $p freq $f");
      }
    }
  }
  
  # check clinsig
  if($v->{clin_sig}) {
    my $valid = get_clinsig_values();
    
    for(split(',', $v->{clin_sig})) {
      ok(exists($valid->{$_}), "clin_sig is valid") or diag("$n clin_sig $_");
    }
  }
  
  # check pubmed
  if($v->{pubmed}) {
    ok($v->{pubmed} =~ /^([0-9]+\,?)+$/, "pubmed looks ok") or diag("$n pubmed ".$v->{pubmed});
  }
}

sub check_rf {
  my $rf = shift;
  my $type = shift;
  check_basics($rf);
  
  ok($rf->isa('Bio::EnsEMBL::Funcgen::'.$type), "RF type check") or diag("expected Bio::EnsEMBL::Funcgen::".$type." got ".$rf);
  
  if($type eq 'MotifFeature') {
    ok($rf->{display_label}, "motif feature has display_label");
    
    if($rf->{binding_matrix}) {
      ok($rf->{binding_matrix}->isa('Bio::EnsEMBL::Funcgen::BindingMatrix'), "binding matrix check");
    }
    
    ok($rf->{_variation_effect_feature_cache} && $rf->{_variation_effect_feature_cache}->{seq}, "motif feature has cached sequence");
  }
}

sub check_tr {
  my $tr = shift;
  check_basics($tr);
  
  ok($tr->isa('Bio::EnsEMBL::Transcript'), "tr type check") or diag("expected Bio::EnsEMBL::Transcript got ".$tr);
  
  ok($tr->{_variation_effect_feature_cache}, "tr has _variation_effect_feature_cache key") or diag(Dumper $tr);
  ok($tr->{_variation_effect_feature_cache}->{mapper} && $tr->{_variation_effect_feature_cache}->{mapper}->isa('Bio::EnsEMBL::TranscriptMapper'), "tr has mapper") or diag(Dumper $tr->{_variation_effect_feature_cache});
  
  # protein coding
  if($tr->{biotype} eq 'protein_coding') {
    $config->{current}->{found_protein_coding} = 1;
    
    ok($tr->{_variation_effect_feature_cache}->{translateable_seq}, "tr has translateable_seq") or diag(Dumper $tr);
    ok($tr->{_variation_effect_feature_cache}->{peptide}, "tr has peptide") or diag(Dumper $tr);
    
    if($tr->{_variation_effect_feature_cache}->{protein_function_predictions}) {
      $config->{current}->{'found_'.$_} = 1 for map {s/\_.+//g; $_} keys %{$tr->{_variation_effect_feature_cache}->{protein_function_predictions}};
    }
  }
}

sub parse_variation {
  my $config = shift;
  my $line = shift;
  
  my @cols = @{$config->{current}->{variation_cols}};
  my @data = split / |\t/, $line;
  
  # assumption fix for old cache files
  if(scalar @data > scalar @cols) {
    push @cols, ('AFR', 'AMR', 'ASN', 'EUR');
  }
  
  my %v = map {$cols[$_] => $data[$_] eq '.' ? undef : $data[$_]} (0..$#data);
  
  $v{failed}  ||= 0;
  $v{somatic} ||= 0;
  $v{end}     ||= $v{start};
  $v{strand}  ||= 1;
  
  return \%v;
}

sub get_clinsig_values {
  if(!exists($config->{current}->{clinsig_values})) {
    my $sth = $config->{current}->{vdb}->prepare(qq{
      DESCRIBE variation
    });
    
    $sth->execute();
    
    my $v = $sth->fetchall_hashref('Field')->{clinical_significance}->{Type};
    $v = (split(/\(|\)/, $v))[1];
    my %valid = map {$_ => 1} map {$_ =~ s/\'//g; $_ =~ s/\s+/\_/g; $_} split(/\,\s*/, $v);
    
    $sth->finish();
    
    $config->{current}->{clinsig_values} = \%valid;
  }
  
  return $config->{current}->{clinsig_values};
}

sub usage {
  
  print qq{
Run this script to check VEP cache installations.

Flags: (defaults follow in parentheses)

--dir | -d           Root directory of VEP caches (\$HOME/.vep/)
--version | -v       Version to check (current API version)
--cache_version | -c Cache version, can be EG version (current API version)
--species | -s       Species to check (all species found on hosts)

--hosts | -h      Database host(s), comma-separated (ens-staging1,ens-staging2)
--user | -u       Database username (ensro)
--password        Database password
--port | -p       Database port (3306)

--no_fasta | -nf  Don't look for and check FASTA file
--random | -r [n] Check content of fraction n cache per chromosome (0)
                  Use --random 1 to check everything!!!
--max_vars | -m   Maximum number of variants to check per cache file (100)

Flags may also be set as env variables prefixed HC_VEP_ e.g.

export HC_VEP_VERSION=78
export HC_VEP_SPECIES=homo_sapiens
};

  return 0;
}
