#!/usr/bin/env perl

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2026] EMBL-European Bioinformatics Institute
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

# Script to load the HGVS CAR lookup files into allele_synonym
# Files loaded contain <variation_id><HGVS><CAid>

use strict;
use warnings;
use Bio::EnsEMBL::Registry;
use Data::Dumper;
use File::Basename;
use POSIX qw/strftime/;

my $pipeline_dir = $ARGV[0];
if (! defined $pipeline_dir) {
  die("Please provide a pipeline_dir");
}

if (! -d $pipeline_dir) {
  die("$pipeline_dir is not a directory: $!");
}

my $car_load_dir = join('/', $pipeline_dir, 'car_load');

if (! -d $car_load_dir) {
  die("$car_load_dir is not a directory: $!");
}

#TODO registry file should be a parameter
my $registry_file = "ensembl.registry.load";

my $registry  = 'Bio::EnsEMBL::Registry';
$registry->load_all($registry_file);

my $species = 'homo_sapiens';
# Get db adaptors to the source and destination databases
my $test_dba = $registry->get_DBAdaptor($species,'variation') or die ("Could not get a DBAdaptor for the variation");

my $dbc = $test_dba->dbc();
$dbc->reconnect_when_lost();

# open a logfile
my $logfile = "car_load.log";
open(my $logh, ">", "$pipeline_dir/$logfile") or die("Unable to open $logfile: $!");
my $run_time = strftime('%Y%m%d_%H%M%S',localtime);
print $logh $run_time, "\n";

# Open the the car_load_dir and get a list of all the JSON fileparse
# hgvs-17-52538080-57315268-load.tab
my @load_files;
print $logh "car_load_dir = $car_load_dir\n";
opendir(LU_DIR, $car_load_dir) or die("Unable to open $car_load_dir: $!");
while (my $file = readdir(LU_DIR))  {
   next if ($file =~ /^\./);
   next if ($file !~ /hgvs-.*-\d{1,}-\d{1,}-load.tab/);
   push @load_files, $file;
 }
close(LU_DIR);

print $logh "Processing ", scalar(@load_files), " load files\n";

for my $file (@load_files) {
  # Check that the load file exists
  print $logh "file: $file\n";
  my $load_file = join("/", $car_load_dir, $file);
  if (! -e $load_file) {
      print $logh "ERROR: load file ($load_file) does not exist\n";
      next;
  }
  if (! -s $load_file) {
      print $logh "ERROR: load file ($load_file) is empty\n";
      next;
  }
  load_car($load_file, $logh, $dbc);
}
print $logh strftime('%Y%m%d_%H%M%S',localtime), "\n";
close($logh);

sub load_car {
  my $load_file = shift;
  my $log_fh = shift;
  my $dbc = shift;

  my $line_count = 0;
  my $rows_loaded = 0;
  my $same = "yes";

  if (! -e $load_file) {
    die("file ($load_file) does not exist");
  }
  my($filename, $dirs, $suffix) = fileparse($load_file, qr/\.[^.]*/);

  $line_count = `wc -l < $load_file`;
  die "wc failed: $?" if $?;
  chomp($line_count);

  my $sql = qq{
    LOAD DATA LOCAL INFILE '$load_file'
    INTO TABLE allele_synonym FIELDS TERMINATED BY '\\t'
    (variation_id, hgvs_genomic, name)
  };
  $rows_loaded = $dbc->do($sql);

  if ($line_count != $rows_loaded) {
    $same = "no";
  }
  print $log_fh join("\t", $filename, $line_count, $rows_loaded, $same), "\n";
}
