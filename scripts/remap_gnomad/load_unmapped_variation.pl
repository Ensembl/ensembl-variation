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
use Bio::EnsEMBL::Registry;

my $dir = $ARGV[0];

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_all($dir . '/remap_gnomad/genomes/ensembl.registry');
my $dbh = $registry->get_DBAdaptor('human', 'variation')->dbc->db_handle;

my $fh = FileHandle->new("$dir/write_gnomad_exomes_unmapped.txt", 'r');

while (<$fh>) {
  chomp;
  my ($seq_region_id, $start, $end, $allele_string, $id) = split/\t/;
  $dbh->do(qq{
    INSERT INTO variation_feature(seq_region_id, seq_region_start, seq_region_end, seq_region_strand, allele_string, variation_name)
    VALUES($seq_region_id, $start, $end, 1, '$allele_string', '$id');
  }) or die $dbh->errstr;
}
$fh->close;
