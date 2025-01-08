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

use FileHandle;

my $grch38_mappings = '/gnomad/Genomes/mapping_results/';
my $grch38_ensembl_mappings = '/remap_gnomad/genomes/unique_mappings/';

my $chrom = $ENV{'LSB_JOBINDEX'};
if ($chrom == 23) {
  $chrom = 'X';
}
if ($chrom == 24) {
  $chrom = 'Y';
}

my $vcf_file = $grch38_mappings . "gnomad.genomes.r2.1.sites.grch38.chr$chrom\_noVEP.vcf";

#print STDERR "gunzip $vcf_file.gz\n";
#run_cmd("gunzip $vcf_file.gz");

open(my $vcf_fh, '>>', $vcf_file) or die "Could not open file '$vcf_file' $!";

my $fh = FileHandle->new("$grch38_ensembl_mappings/$chrom.vcf", 'r');
while (<$fh>) {
  print $vcf_fh $_;
}

$fh->close;
close $vcf_fh;
my $tmpdir = '';
my $cmd = "vcf-sort -t $tmpdir < $vcf_file | bgzip > $vcf_file.gz";

print STDERR "$cmd\n";

run_cmd($cmd);

print STDERR "tabix $vcf_file.gz\n";

run_cmd("tabix $vcf_file.gz");

run_cmd("rm $vcf_file");

sub run_cmd {
  my $cmd = shift;
  if (my $return_value = system($cmd)) {
    $return_value >>= 8;
    die "system($cmd) failed: $return_value";
  }
}

