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

# This script processes the JSON files returned from the 
# CAR (clinical allele registry) pipeline.
# It links them to the file from the HGVS dump that 
# contained the HGVS used to lookup the CA id
# The file produced contain
# <variation_id>\t<hgvs>\t<CA_id>
# The file can be loaded into the allele_synonym table
use strict;
use warnings;
use JSON;
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

my $hgvs_dir = join('/', $pipeline_dir, 'hgvs');
my $car_lu_dir = join('/', $pipeline_dir, 'car_lu');
my $car_load_dir = join('/', $pipeline_dir, 'car_load');

my @dir_list = ($hgvs_dir, $car_lu_dir, $car_load_dir);

for my $dir (@dir_list) {
  if (! -d $dir) {
    die("$dir is not a directory: $!");
  }
}
# open a logfile
my $logfile = "process_lu_car.log";
open(my $logh, ">", "$pipeline_dir/$logfile") or die("Unable to open $logfile: $!");
my $run_time = strftime('%Y%m%d_%H%M%S',localtime);
print $logh $run_time, "\n";

# Open the the car_lu_dir and get a list of all the JSON files
# hgvs-5-146688061-151578662-car-lu.json
my @json_files;
print $logh "car_lu_dir = $car_lu_dir\n";
opendir(LU_DIR, $car_lu_dir) or die("Unable to open $car_lu_dir: $!");
while (my $file = readdir(LU_DIR))  {
   next if ($file =~ /^\./);
   next if ($file !~ /hgvs-.*-\d{1,}-\d{1,}-car-lu.json/);
   push @json_files, $file;
 }
close(LU_DIR);

print $logh "Processing ", scalar(@json_files), " JSON files\n";

for my $file (@json_files) {
  # Check that the hgvs file exists
  # Take the JSON file remove the car-lu
  my $file_base;
  ($file_base = $file)  =~ s/-car-lu.json//;
  print $logh "file: $file\n";
  my $ar_json_file = join("/", $car_lu_dir, $file);
  my $hgvs_file    = join("/", $hgvs_dir, $file_base . ".tab");
  if (! -e $ar_json_file) {
      print $logh "ERROR: JSON file ($ar_json_file) does not exist\n";
      next;
  }
  if (! -e $hgvs_file) {
      print $logh "ERROR: HGVS file ($hgvs_file) does not exist\n";
      next;
  }
  process_ar($ar_json_file, $hgvs_file, $car_load_dir, $logh);
}
print $logh strftime('%Y%m%d_%H%M%S',localtime), "\n";
close($logh);

sub process_ar {
  my $ar_json_file = shift;
  my $hgvs_file = shift;

  my $load_dir = shift;
  my $log_fh = shift;
  if (! -e $ar_json_file) {
    die("file ($ar_json_file) does not exist");
  }
  if (! -e $hgvs_file) {
    die("file ($hgvs_file) does not exist");
  }
  if (! -e $load_dir) {
    die("load dir ($load_dir) is not a directory");
  }

  # Set up the output files
  my($filename, $dirs, $suffix) = fileparse($hgvs_file, qr/\.[^.]*/);

  my $load_file = $filename . "-load.tab";
  my $error_file = $filename . "-error.tab";
  my $missing_file = $filename . "-miss.tab";

  open(my $lh, ">","$load_dir/$load_file") or die("unable to open load_file: $!");
  open(my $eh, ">","$load_dir/$error_file") or die("unable to open error_file: $!");
  open(my $mh, ">","$load_dir/$missing_file") or die("unable to open missing_file: $!");

  # Read in the source data
  my @src_data;
  open(my $sh, '<', $hgvs_file) or die("unable to open $hgvs_file: $!");
  @src_data = <$sh>;
  chomp(@src_data);
  close($sh);

  my $src_count = scalar(@src_data);

  my $json_string;
  {
    open(my $fh, "<", $ar_json_file) or die("unable to open $ar_json_file: $!");
    # read in JSON
    local $/ = undef;
    $json_string = <$fh>;
    close($fh);
  }

  my $data = JSON->new->decode($json_string) or throw("ERROR: Failed to parse json $json_string");
  my $json_count = scalar(@$data);

  if ($src_count != $json_count) {
    print $log_fh "ERROR: src_count ($src_count) != json_count ($json_count)\n";
    return;
  }

  my $car_count = 0;
  my $miss_count = 0;
  my $err_count = 0;
  my $len_count = 0;

  for (my $i = 0; $i <= $#src_data; $i++) {
    # need to get the fifth element of the inputLine
    my $hgvs = (split("\t", $src_data[$i]))[5];
    if (! $hgvs) {
      print $log_fh "ERROR: $filename hgvs missing ", (split("\t",$src_data[$i]))[0], "\n";
      $len_count++;
      next;
    }
    if (length($hgvs) > 600) {
      print $log_fh "ERROR: $filename hgvs too long ", (split("\t",$src_data[$i]))[0,5], "\n";
      $len_count++;
      next;
    }
    if (exists $data->[$i]->{'@id'}) {
      my ($car_id) = ($data->[$i]->{'@id'} =~/(CA\d{1,})$/);
      if ($car_id) {
        print $lh join("\t", (split("\t",$src_data[$i]))[0,5], $car_id), "\n";
        $car_count++;
      } else {
        print $mh join("\t", (split("\t",$src_data[$i]))[0,5]), "\n";
        $miss_count++;
      }
    } else {
      $err_count++;
      my $ao = $data->[$i];
      next if (! $ao->{'errorType'});
      print $eh join("\t",
        $src_data[$i],
        $ao->{'inputLine'},
        $ao->{'hgvs'} || '',
        $ao->{'errorType'},
        $ao->{'message'},
        $ao->{'description'}),
        "\n";
    }
  }
  my $total = scalar(@$data);
  print $log_fh join("\t", $filename, $total, $car_count, $miss_count, $err_count, $len_count), "\n";
}
