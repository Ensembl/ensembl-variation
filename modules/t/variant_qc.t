#!/usr/bin/env perl
# Copyright 2013 Ensembl
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


## initial tests for hive QC pipeline

use strict;
use warnings;


use Test::More;
use Data::Dumper;
use FindBin qw($Bin);

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::Utils::Sequence qw(revcomp_tandem);
use Bio::EnsEMBL::Variation::Utils::QCUtils qw( check_four_bases check_illegal_characters remove_ambiguous_alleles find_ambiguous_alleles check_variant_size);

use_ok('Bio::EnsEMBL::Variation::Pipeline::VariantQC::VariantQC_conf');
use_ok('Bio::EnsEMBL::Variation::Pipeline::VariantQC::InitVariantQC');
use_ok('Bio::EnsEMBL::Variation::Pipeline::VariantQC::CheckdbSNPImport');
use_ok('Bio::EnsEMBL::Variation::Pipeline::VariantQC::VariantQC');
use_ok('Bio::EnsEMBL::Variation::Pipeline::VariantQC::UpdatePopulationGenotype');
use_ok('Bio::EnsEMBL::Variation::Pipeline::VariantQC::FlipPopulationGenotype');
use_ok('Bio::EnsEMBL::Variation::Pipeline::VariantQC::FinishVariantQC');
use_ok('Bio::EnsEMBL::Variation::Pipeline::VariantQC::RegisterDBSNPImport');




## tiny test database with plenty of fails needed...


ok( revcomp_tandem("(GT)17/(GT)19") eq "(AC)17/(AC)19", "Utils::Sequence revcomp_tandem" );

ok( check_four_bases("A/T/G/G/C")   eq 1 ,              "Utils::QCUtils check_four_bases positive");
ok( check_four_bases("A/T/-/AA/C")  eq 0 ,              "Utils::QCUtils check_four_bases negative");


my $ambiguous_alleles_pos = find_ambiguous_alleles("M/T" );
ok ($ambiguous_alleles_pos->[0]        eq "M",             "Utils::QCUtils find_ambiguous_alleles positive");

my $ambiguous_alleles_neg = find_ambiguous_alleles("A/T");
ok (! defined $ambiguous_alleles_neg->[0],                "Utils::QCUtils find_ambiguous_alleles negative");


my  $illegal_characters_pos = check_illegal_characters("FISH");
ok( $illegal_characters_pos->[0] eq "FISH",          "Utils::QCUtils illegal_characters positive");

my  $illegal_characters_neg = check_illegal_characters("A/T");
ok( ! defined $illegal_characters_neg->[0] ,          "Utils::QCUtils illegal_characters negative");

my  $illegal_characters_neg2 = check_illegal_characters("M/T");
ok( ! defined $illegal_characters_neg2->[0] ,          "Utils::QCUtils illegal_characters negative (ambig)");


ok( remove_ambiguous_alleles("M/T") eq "[AC]/T",         "Utils::QCUtils remove_ambiguous_alleles positive ");
ok( remove_ambiguous_alleles("A/T") eq "A/T",            "Utils::QCUtils remove_ambiguous_alleles negative ");

ok( check_variant_size(10,12,"A")   eq 0,                "Utils::QCUtils variant size error1");
ok( check_variant_size(10,10,"AA")  eq 0,                "Utils::QCUtils variant size error2");
ok( check_variant_size(10,12,"-")   eq 0,                "Utils::QCUtils variant size error3");
ok( check_variant_size(10,10,"-")   eq 1,                "Utils::QCUtils variant size insertion");
ok( check_variant_size(10,12,"AAA") eq 1,                "Utils::QCUtils variant size del");
ok( check_variant_size(10,10,"A")   eq 1,                "Utils::QCUtils variant size error");


done_testing();

