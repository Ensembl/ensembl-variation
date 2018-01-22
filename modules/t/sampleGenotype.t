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

use Test::Exception;
use Test::More;
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Variation::Individual;
use Bio::EnsEMBL::Variation::Sample;
use Bio::EnsEMBL::Variation::SampleGenotype;
use Bio::EnsEMBL::Variation::Variation;
our $verbose = 0;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');
my $vdba = $multi->get_DBAdaptor('variation');
my $sga = $vdba->get_SampleGenotypeAdaptor;

# test constructor
my $ind = Bio::EnsEMBL::Variation::Individual->new(
  -name => 'test individual',
  -description => 'This is a test individual',
  -gender => 'Male');

my $sample = Bio::EnsEMBL::Variation::Sample->new(
  -name => 'test sample',
  -individual => $ind,
);

my $var = Bio::EnsEMBL::Variation::Variation->new(
  -name => 'rs123',
  -synonyms => {'dbSNP' => ['ss12', 'ss144']},
  -source => 'dbSNP');

my $genotype = ['A','C'];
my $subsnp   ='ss12';

throws_ok { 
  my $sample_gtype = Bio::EnsEMBL::Variation::SampleGenotype->new(
    -genotype => $genotype,
    -variation => 'variation',
    -sample => $sample,
    -subsnp => $subsnp  
  );
} qr/Bio::EnsEMBL::Variation::Variation argument expected/, 'Throw on wrong argument for new';

throws_ok { 
  my $sample_gtype = Bio::EnsEMBL::Variation::SampleGenotype->new(
    -genotype => $genotype,
    -variation => $var,
    -sample => 'sample',
    -subsnp => $subsnp  
  );
} qr/Bio::EnsEMBL::Variation::Sample argument expected/, 'Throw on wrong argument for new';

my $sample_gtype = Bio::EnsEMBL::Variation::SampleGenotype->new(
  -genotype => $genotype,
  -variation => $var,
  -sample => $sample,
  -subsnp => $subsnp  
);

ok($sample_gtype->sample()->name() eq $sample->name(), "sample name");
ok($sample_gtype->variation()->name()  eq $var->name(), "var name");

ok($sample_gtype->allele1() eq 'A', "allele 1");
ok($sample_gtype->allele2() eq 'C', "allele 2");

ok($sample_gtype->genotype() eq $genotype, "genotype");
my $string = join '|', @{$genotype};
ok($sample_gtype->genotype_string() eq $string, "genotype string");

ok($sample_gtype->subsnp() eq $subsnp, "subsnp");

ok($sample_gtype->ambiguity_code() eq 'M', "ambiguity code"); 

# test getter/setters
ok(test_getter_setter($sample_gtype, 'allele1', 'TT'), "get/set allele 1");
ok(test_getter_setter($sample_gtype, 'allele2', '-'),  "get/set allele 2");

my $ind2 = Bio::EnsEMBL::Variation::Individual->new(
  -name => 'test individual 2',
  -description => 'This is a second test individual',
  -gender => 'Female'
);

my $sample2 = Bio::EnsEMBL::Variation::Sample->new(
  -name => 'test sample 2',
  -description => 'This is a second test sample',
  -individual => $ind2,
);

my $var2 = Bio::EnsEMBL::Variation::Variation->new(
  -name => 'rs332',
  -synonyms => {'dbSNP' => ['ss44', 'ss890']},
  -source => 'dbSNP'
);

ok(test_getter_setter($sample_gtype, 'sample', $sample2), "get/set sample");
ok(test_getter_setter($sample_gtype, 'variation', $var2),  "get/set variation");

throws_ok { $sample_gtype->sample('sample'); } qr/Bio::EnsEMBL::Variation::Sample argument expected/, 'Throw on wrong argument for sample';

$sample_gtype = Bio::EnsEMBL::Variation::SampleGenotype->new(
  -genotype => $genotype,
  -variation => $var,
  -subsnp => $subsnp,  
  -adaptor => $sga
);
$sample_gtype->{sample_id} = 2105; 
$sample = $sample_gtype->sample();
ok($sample->name eq 'NA18635', 'get sample');
done_testing();
