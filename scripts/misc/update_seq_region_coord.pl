#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2023] EMBL-European Bioinformatics Institute
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

use strict;
use warnings;

use Bio::EnsEMBL::DBSQL::DBConnection;
use Getopt::Long;
use Bio::EnsEMBL::Utils::CliHelper;

my $script_opts = [{args => ['host', 'dbhost', 'h'], type => '=s'},
        {args => ['port', 'dbport', 'P'], type => ':i'},
        {args => ['user', 'dbuser', 'u'], type => '=s'},
        {args => ['pass', 'dbpass', 'p'], type => ':s'},
        {args => ['dbname',  'D'],         type => ':s'},
        {args => ['version'],    type => ':i'},
        {args => ['dry_run'], type=>'!'}, ],
        {args => ['truncate'], type=>'!'}, ],;

if (scalar(@ARGV) == 0) {
  usage();
  exit 0;
}

my $cli_helper = Bio::EnsEMBL::Utils::CliHelper->new();
my $optsd = get_opts();

my $opts = $cli_helper->process_args( $optsd, \&usage );

my ($variation_dbas, $core_dbas);
if (defined $opts->{dbname}) {
  die "Variation database required" unless ($opts->{dbname} =~ m/variation/);
  $variation_dbas = get_dbas($opts);
} else {
  my $version = $opts->{version};
  die "Version argument is required" unless ($version);
  $opts->{pattern} = 'variation_'.$version.'_\d+$';
  $variation_dbas = get_dbas($opts);
}

foreach my $variation_dbname (keys %$variation_dbas) {
  print "variation db: $variation_dbname\n";
  my $core_dbname = $variation_dbname;
  $core_dbname =~ s/variation/core/;
  print "core db: $core_dbname\n";
  my $vdba = $variation_dbas->{$variation_dbname};
  next unless ($vdba);
  
  # Coord system need to be truncated?
  if ($opts->{truncate}) {
    my $sql_truncate = qq{
      TRUNCATE table coord_system
    };
    my $trunc = $vdba->dbc()->sql_helper()->execute(
                                         -SQL => $sql_truncate);
    print "Total of $trunc rows in coord_system set\n";
  }

  # Are there foreign key failures between seq_region and coord_system?
  my $fk_count = count_fk_failures($vdba);
  if ($fk_count == 0) {
    print "$variation_dbname: No foreign key errors between seq_region and coord_system tables\n";
    print "$variation_dbname: No update needed to seq_region table\n";
    next;
  } else {
    print "$variation_dbname: $fk_count foreign key errors between seq_region and coord_system tables\n";
  }
  
  # Update the seq_region
  print "$variation_dbname: START seq_region update\n";
  my $sql_update = qq{
      UPDATE seq_region vsr, 
         $core_dbname.seq_region csr
      SET vsr.coord_system_id = csr.coord_system_id
      WHERE vsr.seq_region_id = csr.seq_region_id
    };
  
  if ($opts->{dry_run}) {
    print "$variation_dbname: $sql_update\n";
  } else {
    my $num_update = $vdba->dbc()->sql_helper()->execute_update(
                                         -SQL => $sql_update);
    if ($num_update) {
      print "$variation_dbname: $num_update seq_region records updated\n";
    } else {
      print "$variation_dbname: no seq_region records updated\n";
    }
  }
  print "$variation_dbname: END seq_region update\n";
  
  # Update the coord_system
  print "$variation_dbname: START coord_system update\n";
  $sql_update = qq{
    INSERT IGNORE INTO coord_system
    (`coord_system_id`, `species_id`, `name`, `version`, `rank`, `attrib`)
    SELECT coord_system_id, species_id, name, version, rank, attrib
    FROM $core_dbname.coord_system cs
    WHERE coord_system_id IN
      (SELECT DISTINCT coord_system_id 
       FROM seq_region
       );
  };
  
  if ($opts->{dry_run}) {
    print "$variation_dbname: $sql_update\n";
  } else {
    my $num_update = $vdba->dbc()->sql_helper()->execute_update(
                                         -SQL => $sql_update);
    print "$variation_dbname: $num_update coord_system records added\n";
  }
  
  print "$variation_dbname: END coord_system update\n";
    
  # Are there foreign key failures remaining?
  $fk_count = count_fk_failures($vdba);
  if ($fk_count == 0) {
      print "$variation_dbname: No foreign key errors between seq_region and coord_system tables\n";
  } else {
    print "$variation_dbname: $fk_count foreign key errors between seq_region and coord_system tables\n";
  }
}

sub get_dbas { 
  my $opts = shift;
  my $dbas = {};
  # use the command line options to get an array of database details
  for my $db_args ( @{ $cli_helper->get_dba_args_for_opts( $opts, 0 ) } ) {
    # use the args to create a DBA
    my $dba = new Bio::EnsEMBL::DBSQL::DBAdaptor(%$db_args);
    my $dbname = $dba->dbc()->dbname();
    $dbas->{$dbname} = $dba;
  }
  return $dbas;
}

sub get_opts {
  my @dba_opts = map {
    my $opt = join '|', @{$_->{args}};
    $opt . $_->{type};
  } @{$script_opts};
  return \@dba_opts;
}

sub usage {
  print <<USAGE_END;
USAGE:
  $0 --dbhost=XhostX [--dbport=3306] \\
  \t--dbuser=XuserX --dbpass=XXX \\
  $0 \$(XserverX details script_db) --dbname XdbnameX
  $0 \$(XserverX details script_db) --dbname XdbnameX --dry_run
  $0 \$(XserverX details script_db) --version Xrelease-versionX
  $0 --help
  --dry_run     Print update statements
  --help        Displays this help text.

This script will update 
- the variation seq_region.coord_system_id with core seq_region.coord_system_id
- the variation coord_system table with coord_system records
  that are in variation seq_region.coord_system_id from core coord_system

The script expects the core and variation database to be on the same db host.

USAGE_END

}

sub count_fk_failures {
  my ($vdba) = @_;
  my $helper = $vdba->dbc()->sql_helper();
  
  my $sql = qq{
     SELECT COUNT(*) 
     FROM
       seq_region t1 LEFT JOIN coord_system t2 
        ON t1.coord_system_id = t2.coord_system_id
     WHERE t1.coord_system_id IS NOT NULL AND t2.coord_system_id IS NULL
   };
    
   my $count = $helper->execute_single_result(-SQL => $sql);
   return $count;
}
