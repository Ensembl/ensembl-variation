#!/usr/bin/env perl

=head1 LICENSE

Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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
use Bio::EnsEMBL::Registry;
BEGIN {
  $| = 1;
	use Test::More;
}

my $config = {};

GetOptions(
  $config,
  'help|h',                  # displays help message
  
  'dir|d=s',                 
  'version|v=i',
  'species|s=s',
  
  'host|h=s',                  # DB options
  'user|u=s',
  'password|p=s',
  'port|P=i',
) or die "ERROR: Failed to parse command-line flags\n";

do { usage(); exit(0); } if $config->{help};

$config->{reg} = 'Bio::EnsEMBL::Registry';

$config->{dir}     ||= $ENV{HOME}.'/.vep/';
$config->{version} ||= $config->{reg}->software_version;
$config->{species} ||= 'all';

$config->{host} ||= 'ens-staging1,ens-staging2';
$config->{user} ||= 'ensro';
$config->{port} ||= 3306;

$config->{cache_region_size} = 1e6;

my %match_species = ();
if($config->{species} ne 'all') {
  $match_species{$_} = 1 for split(/\,/, $config->{species});
}

foreach my $host(split /\,/, $config->{host}) {
  
  my $species_list = get_species_list($config, $host);
  
  if($config->{species} eq 'all') {
    $match_species{$_} = 1 for @$species_list;
  }
  
  $config->{reg}->load_registry_from_db(
    -host       => $host,
    -user       => $config->{user},
    -pass       => $config->{password},
    -port       => $config->{port},
    -db_version => $config->{version},
  );
  
  foreach my $species(@$species_list) {
    next unless $match_species{$species};
    
    print "\nChecking $species\n";
    
    my $sp = $species;
    my $refseq = $sp =~ s/\_refseq//;
    
    my $has_var = has_variation($config, $sp);
    my $has_reg = has_regulation($config, $sp);
    
    # check directory structure
    print " - checking directory structure\n";
    
    my $dir = $config->{dir}.'/'.$species.'/'.$config->{version};
    
    ok(-d $config->{dir}.'/'.$species, "\[$species\] species dir exists");
    ok(-d $dir, "\[$species\] version dir exists");
    
    
    # check info file
    print " - checking info.txt\n";
    
    ok(-e $dir.'/info.txt', "\[$species\] info.txt exists");
    ok(open(IN, $dir.'/info.txt'), "\[$species\] info.txt readable");
    
    my ($v, $r, $tabix);
    
    while(<IN>) {
      ok(/all/, "--build all used") if /^build/;
      $v = $_ if /^variation_cols/;
      $r = $_ if /^regulatory/;
      $tabix = $_ if /^var_type/;
    }
    close IN;
    
    ok(defined($v), "\[$species\] variation_cols defined in info.txt") if $has_var;
    ok(defined($r), "\[$species\] regulation defined in info.txt") if $has_reg;
    
    
    # check chromosome dirs exist
    print " - checking chromosome dirs";
    
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
    $counts{$sr_id} = $count while $sth->fetch;
    $sth->finish;
    
    foreach my $slice(@slices) {
      #printf("\r - checking chromosome dirs %i / %i", ++$i, scalar @slices);
      
      my $chr = $slice->seq_region_name;
      
      next unless $counts{$slice->get_seq_region_id};
      
      push @missing_dirs, $chr unless -d $dir.'/'.$chr;
      
      # initial region
      my $start = 1 + ($config->{cache_region_size} * int($slice->start / $config->{cache_region_size}));
      my $end   = ($start - 1) + $config->{cache_region_size};
      
      while($start < $slice->end) {
        
        my @types = ('');
        push @types, '_var' if $has_var && !$tabix;
        push @types, '_reg' if $has_reg;
        
        foreach my $type(@types) {
          my $file = $dir.'/'.$chr.'/'.$start.'-'.$end.$type.'.gz';
          push @missing_files, $file unless -e $file;
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
        }
      }
    }
    print "\n";
    
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
    
    
    # check FASTA file
    print " - checking FASTA file\n";
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

print "\n\nDone testing\n";
done_testing();

sub get_species_list {
  my $config = shift;
  my $host   = shift;

  my $connection_string = sprintf(
    "DBI:mysql(RaiseError=>1):host=%s;port=%s;db=mysql",
    $host,
    $config->{port}
  );
  
  # connect to DB
  $config->{dbc} = DBI->connect(
    $connection_string, $config->{user}, $config->{password}
  );
  
  my $version = $config->{version};

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

  foreach my $current_db_name (@dbs) {
    
    # special case otherfeatures
    # check it has refseq transcripts
    if($current_db_name =~ /otherfeatures/) {
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
        push @species, $current_db_name;
      }
    }
    
    else {
      $sth = $config->{dbc}->prepare("select meta_value from ".$current_db_name.".meta where meta_key='species.production_name';");
      $sth->execute();
      my $current_species = $sth->fetchall_arrayref();
      
      my @flattened_species_list = sort map { $_->[0] } @$current_species;
      
      if(@flattened_species_list) {
        push @species, @flattened_species_list;
      }
      else {
        $current_db_name =~ s/^([a-z]+\_[a-z,1-9]+)(\_[a-z]+)?(.+)/$1$2/;
        $current_db_name =~ s/\_core$//;
        push @species, $current_db_name;
      }
    }
  }

  return \@species;
}

sub has_variation {
  my ($config, $species) = @_;
  
  my $version = $config->{version};
  
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
  
  my $version = $config->{version};
  
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
    SELECT meta_value
    FROM $db.meta
    WHERE meta_key = 'regbuild.version'
  })->[0]->[0];
  
  return defined($v) && $v > 0;
}

sub usage {
  
  print qq{
Run this script to check VEP cache installations.

Flags: (defaults follow in parentheses)

--dir | -d        Root directory of VEP caches (\$HOME/.vep/)
--version | -v    Version to check (current API version)
--species | -s    Species to check (all species found on hosts)

--hosts | -h      Database host(s), comma-separated (ens-staging1,ens-staging2)
--user | -u       Database username (ensro)
--password | -p   Database password
--port | -P       Database port (3306)
};

  return 0;
}
