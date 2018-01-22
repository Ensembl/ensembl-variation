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

foreach my $cdba (@$cdbas) {
  my $dbname = $cdba->dbc->dbname;
  my $species = $cdba->species;
  next if (!$vdbas->{$species});
  my $dbh = $cdba->dbc->db_handle;
  my $max_mapping_set_id = get_max_mapping_set_id($dbh, $dbname); 
  my $id_mapping = {};
  my $sth = $dbh->prepare("SELECT external_seq_region_id, internal_seq_region_id FROM seq_region_mapping WHERE mapping_set_id=$max_mapping_set_id;");
  $sth->execute();
  while (my @row = $sth->fetchrow_array) {
    my $external_seq_region_id = $row[0]; # previous release
    my $internal_seq_region_id = $row[1]; # current release
    $id_mapping->{$external_seq_region_id} = $internal_seq_region_id;     
  }
  $sth->finish();
  
  my $vdba = $vdbas->{$species};
  my $vdbh = $vdba->dbc->db_handle;
  foreach my $prev_seq_region_id ( keys %$id_mapping) {
    my $new_seq_region_id = $id_mapping->{$prev_seq_region_id};
    if ($config->{dry_run}) {
      print "For $dbname: Update seq_region SET seq_region_id=$new_seq_region_id WHERE seq_region_id=$prev_seq_region_id\n";
    } else {
      $vdbh->do("Update seq_region SET seq_region_id=$new_seq_region_id WHERE seq_region_id=$prev_seq_region_id") or die $dbh->errstr;
    }
  }
}

sub get_max_mapping_set_id {
  my ($dbh, $dbname) = @_;
  my $sth_mapping = $dbh->prepare("select max(mapping_set_id) from $dbname.mapping_set");
  $sth_mapping->execute();
  my ($max_mapping_set_id) = $sth_mapping->fetchrow_array();
  if (!defined $max_mapping_set_id) { return 0; }
  return $max_mapping_set_id;
}

sub usage {
  print qq{
  Usage: perl update_seq_region_ids.pl -registry [registry_file] [OPTIONS]
  Update seq_region_ids between releases. Check if any seq_region_ids have changed since the last release and update them if needed. 
  Options:
    -help    Print this message
    -dry_run Print update statements
  } . "\n";
  exit(0);
}
