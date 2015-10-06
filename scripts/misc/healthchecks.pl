#!/usr/bin/env perl

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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
  'species=s',
  'help!',
) or die "Error: Failed to parse command line arguments\n";

usage() if ($config->{help});

die ('A registry file is required (--registry)') unless (defined($config->{registry}));

my $registry = 'Bio::EnsEMBL::Registry';   
$registry->load_all($config->{registry});
$config->{species} ||= 'human';
my $species = $config->{species};
my $vdba = $registry->get_DBAdaptor($species, 'variation');
my $dbh = $vdba->dbc->db_handle;

my $arguments = [
  ["variation", "variation_id", "somatic", "variation_feature", "variation_id", "somatic"],
  ["variation_feature", "variation_feature_id", "somatic", "transcript_variation", "variation_feature_id", "somatic"],
  ["variation", "variation_id", "display", "variation_feature", "variation_id", "display"],
  ["variation_feature", "variation_feature_id", "display", "transcript_variation", "variation_feature_id", "display"],
];

foreach my $argument_list (@$arguments) {
  check_for_bad_denormalization(@$argument_list);
}

sub check_for_bad_denormalization {
  my ($table1, $col1, $col1d, $table2, $col2, $col2d) = @_;

  my $sql = "$table1, $table2 WHERE $table1.$col1 = $table2.$col2 AND $table1.$col1d != $table2.$col2d";
  my $useful_sql = "SELECT $table1.$col1, $table1.$col1d, $table2.$col2, $table2.$col2d FROM $sql"; 

  print STDERR "Running SELECT count(*) FROM $sql\n";

  my $sth = $dbh->prepare("SELECT count(*) FROM $sql", {mysql_use_result => 1});
  $sth->execute();
  while (my @row = $sth->fetchrow_array) {
    my $count = $row[0];
    if ($count > 0) {
      print STDERR "FAILED $table1 -> $table2 on the denormalization of $col1d using the FK $col1\n";
      print STDERR "FAILURE DETAILS: $count $col1d entries are different in $table1 and $table2\n";
      print STDERR "USEFUL SQL: $useful_sql\n";
    }
  }
  $sth->finish();
}

sub usage {
  print qq{
  Usage: perl healthchecks.pl -registry [registry_file] [OPTIONS]
  Run long running healthchecks separately from overnight healthchecks.

  Options:
    -species     Default is human
    -help        Print this message
  } . "\n";
  exit(0);
}
