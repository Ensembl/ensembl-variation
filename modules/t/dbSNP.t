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


use Test::More;

use Data::Dumper;

use FindBin qw($Bin);

use Bio::EnsEMBL::Variation::Utils::dbSNP qw(get_alleles_from_pattern);

my @test_patterns = ( ["(AGAC)25/26/27/28/29/30/31/32/33/34/35/36/37/38/39/40/41/42/43/44/", 
                      "(AGAC)25/(AGAC)26/(AGAC)27/(AGAC)28/(AGAC)29/(AGAC)30/(AGAC)31/(AGAC)32/(AGAC)33/(AGAC)34/(AGAC)35/(AGAC)36/(AGAC)37/(AGAC)38/(AGAC)39/(AGAC)40/(AGAC)41/(AGAC)42/(AGAC)43/(AGAC)44"],
                      ["(AT)7/8", "(AT)7/(AT)8"],
		      ["T/C", "T/C"], 
		      ["-/A", "-/A"],
		      ["(FB66D06.X1)",  "(FB66D06.X1)"],
                      ["(TG)20/21/A/G", "(TG)20/(TG)21/A/G"],
                      ["(TTA)1/6/7/-",  "TTA/(TTA)6/(TTA)7/-"]
		      );


foreach my $t (@test_patterns){

  my $als = get_alleles_from_pattern($t->[0]);

  my $ens_string = join "/", @{$als}; 
#  warn "\n$ens_string\n";

  ok( $t->[1] eq $ens_string , "[ testing: $t->[0] -> $t->[1] ]");


}

done_testing();
