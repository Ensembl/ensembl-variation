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

use Test::More;
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Variation::Individual;
use Bio::EnsEMBL::Variation::IndividualGenotype;
use Bio::EnsEMBL::Variation::Variation;
our $verbose = 0;

my $ind = Bio::EnsEMBL::Variation::Individual->new(
  -name => 'test individual',
  -description => 'This is a test individual',
  -gender => 'Male'
);

my $var = Bio::EnsEMBL::Variation::Variation->new(
  -name => 'rs123',
  -synonyms => {'dbSNP' => ['ss12', 'ss144']},
  -source => 'dbSNP'
);

my $genotype = ['A','C'];
my $subsnp   ='ss12';

my $ind_gtype = Bio::EnsEMBL::Variation::IndividualGenotype->new(
  -genotype => $genotype,
  -variation => $var,
  -individual => $ind,
  -subsnp     => $subsnp  
);

ok($ind_gtype->individual()->name() eq $ind->name(), "ind name");

done_testing();
