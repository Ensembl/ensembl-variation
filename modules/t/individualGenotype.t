# Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Variation::Individual;
use Bio::EnsEMBL::Variation::IndividualGenotype;
use Bio::EnsEMBL::Variation::Variation;
our $verbose = 0;

# test constructor
my $ind = Bio::EnsEMBL::Variation::Individual->new(
  -name => 'test individual',
  -description => 'This is a test individual',
  -gender => 'Male');

my $var = Bio::EnsEMBL::Variation::Variation->new(
  -name => 'rs123',
  -synonyms => {'dbSNP' => ['ss12', 'ss144']},
  -source => 'dbSNP');

my $genotype = ['A','C'];
my $subsnp   ='ss12';

my $ind_gtype = Bio::EnsEMBL::Variation::IndividualGenotype->new(
  -genotype => $genotype,
  -variation => $var,
  -individual => $ind,
  -subsnp     => $subsnp  
);

ok($ind_gtype->individual()->name() eq $ind->name(), "ind name");
ok($ind_gtype->variation()->name()  eq $var->name(), "var name");

ok($ind_gtype->allele1() eq 'A', "allele 1");
ok($ind_gtype->allele2() eq 'C', "allele 2");

ok($ind_gtype->genotype() eq $genotype, "genotype");
my $string = join '|', @{$genotype};
ok($ind_gtype->genotype_string() eq $string, "genotype string");

ok($ind_gtype->subsnp() eq $subsnp, "subsnp");

ok($ind_gtype->ambiguity_code() eq 'M', "ambiguity code"); 

# test getter/setters
ok(test_getter_setter($ind_gtype, 'allele1', 'TT'), "get/set allele 1");
ok(test_getter_setter($ind_gtype, 'allele2', '-'),  "get/set allele 2");

my $ind2 = Bio::EnsEMBL::Variation::Individual->new
  (-name => 'test individual 2',
   -description => 'This is a second test individual',
   -gender => 'Female');

my $var2 = Bio::EnsEMBL::Variation::Variation->new
  (-name => 'rs332',
   -synonyms => {'dbSNP' => ['ss44', 'ss890']},
   -source => 'dbSNP');

ok(test_getter_setter($ind_gtype, 'individual', $ind2), "get/set individual");
ok(test_getter_setter($ind_gtype, 'variation', $var2),  "get/set variation");


done_testing();
