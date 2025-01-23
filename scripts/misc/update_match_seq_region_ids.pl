#!/usr/bin/env perl

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2025] EMBL-European Bioinformatics Institute
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

use Bio::EnsEMBL::Registry;

use DBI qw(:sql_types);
use Getopt::Long;

usage() if (!scalar(@ARGV));

my $config = {};

GetOptions(
  $config,
  'registry=s',
  'dry_run',
  'help!',
) or die "Error: Failed to parse command line arguments\n";

usage() if ($config->{help});

die ('A registry file is required (--registry)') unless (defined($config->{registry}));

my $registry = 'Bio::EnsEMBL::Registry';   
$registry->load_all($config->{registry});

my $cdbas = $registry->get_all_DBAdaptors(-group => 'core');
my $vdbas_tmp = $registry->get_all_DBAdaptors(-group => 'variation');
my $vdbas = {};
foreach my $vdba (@$vdbas_tmp) {
  my $species = $vdba->species;
  $vdbas->{$species} = $vdba;
}

my @tables = (
  "compressed_genotype_region",
  "phenotype_feature",
  "read_coverage",
  "variation_feature",
  "structural_variation_feature"
  );

foreach my $cdba (@$cdbas) {
  my $species = $cdba->species;
  die ("Species not found") if (!$vdbas->{$species});
  my $dbh = $cdba->dbc->db_handle;
  my $cd_dbname = $cdba->dbc->dbname;
  my $id_mapping = {};
  my $sth = $dbh->prepare("SELECT seq_region_id, name, coord_system_id FROM seq_region;");
  $sth->execute();
  while (my @row = $sth->fetchrow_array) {
    my $external_seq_region_id = $row[0];
    my $external_seq_region_name = $row[1] . "|" . $row[2];
    $id_mapping->{$external_seq_region_name} = $external_seq_region_id;    
  }
  $sth->finish();
  
  my $vdba = $vdbas->{$species};
  my $vdbh = $vdba->dbc->db_handle;

  my $vd_dbname = $vdba->dbc->dbname;
  my %vd_mapping = ();
  $sth = $vdbh->prepare("SELECT seq_region_id, name, coord_system_id FROM seq_region;");
  $sth->execute();
  while (my @row = $sth->fetchrow_array) {
    my $internal_seq_region_id = $row[0];
    my $internal_seq_region_name = $row[1] . "|" . $row[2];

    # Check if Variation and Core are equal in a name
    die "ERROR: '$internal_seq_region_name' is not listed in $cd_dbname.seq_region name column\n." unless defined($id_mapping->{$internal_seq_region_name});

    $vd_mapping{$internal_seq_region_name} = $internal_seq_region_id;
  }
  $sth->finish();

  # Core vs. Variation seq_region_id's
  die ("Variation vs Core: all seq_region ids match\n") unless (%vd_mapping);

  # Remove old seq_region_id from Variation db and create a new one based on core
  unless (defined($config->{dry_run})) {
    $vdbh->do("ALTER TABLE seq_region drop seq_region_id;") or die $dbh->errstr;
    $vdbh->do("ALTER TABLE seq_region ADD seq_region_id INT UNSIGNED NOT NULL FIRST") or die $dbh->errstr;
  }

  print "Checking if table exists and is populated in Variation ... \n";
  # Table check if is populated and exists
  my $need_ids = {};
  foreach my $table (@tables){
    my $count = check_if_exists($vdba, $table);
    if ($count == 0) {
      print "Skipping empty table '$table'\n";
      next;
    }

    # Get unique ids from table
    my @ids = get_ids($vdbh, $table);
    $need_ids->{$table} = \@ids;
  }
  print "OK\n";

  foreach my $prev_seq_region ( keys %$id_mapping) {
    my $new_seq_region_id = $id_mapping->{$prev_seq_region};
    my $old_seq_region_id = $vd_mapping{$prev_seq_region};
    my ($prev_seq_region_name, $prev_seq_region_coord);
    my @psr = split(/\|/, $prev_seq_region);
    $prev_seq_region_name = $psr[0];
    $prev_seq_region_coord = $psr[1];

    # Skip if old and new are the same
    next if (!$old_seq_region_id);

    if ($config->{dry_run}) {
      print "Update seq_region SET seq_region_id=$new_seq_region_id WHERE name='$prev_seq_region_name' AND coord_system_id=$prev_seq_region_coord\n";
      foreach my $table (@tables) {
        next if ( $new_seq_region_id eq $old_seq_region_id);
        next if ( ! grep $_ eq $old_seq_region_id, @{$need_ids->{$table}});
        print "Update $table SET seq_region_id=$new_seq_region_id WHERE seq_region_id=$old_seq_region_id\n";
      }
    } else {
      $vdbh->do("Update seq_region SET seq_region_id=$new_seq_region_id WHERE name='$prev_seq_region_name' AND coord_system_id=$prev_seq_region_coord") or die $dbh->errstr;
      foreach my $table (@tables) {
        next if ( $new_seq_region_id eq $old_seq_region_id);
        next if ( ! grep $_ eq $old_seq_region_id, @{$need_ids->{$table}});
        $vdbh->do("Update $table SET seq_region_id=$new_seq_region_id WHERE seq_region_id=$old_seq_region_id") or die $dbh->errstr;
      }
    }
  }

  unless (defined($config->{dry_run})) {
    $vdbh->do("ALTER TABLE seq_region ADD PRIMARY KEY (seq_region_id);") or die $dbh->errstr;
  }

}

sub usage {
  print qq{
  Usage: perl update_seq_region_ids.pl -registry [registry_file] [OPTIONS]
  Update seq_region_ids between core and variation.
  Options:
    -help    Print this message
    -dry_run Print update statements
  } . "\n";
  exit(0);
}

sub check_if_exists {
  my ($vdba, $table) = @_;
  my $helper = $vdba->dbc()->sql_helper();
  
  my $sql = qq{
     SELECT COUNT(*) 
     FROM
       $table
   };
    
   my $count = $helper->execute_single_result(-SQL => $sql);
   return $count;
}

sub get_ids {
  my ($vdbh, $table) = @_;
  
  my @ids = ();
  my $sth = $vdbh->prepare("SELECT DISTINCT seq_region_id FROM $table;");
  $sth->execute();
  while (my @row = $sth->fetchrow_array) {
    push(@ids, $row[0]);
  }
  $sth->finish();
  
  return @ids;
}