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
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Variation::Individual;
use Bio::EnsEMBL::Variation::Sample;
use Bio::EnsEMBL::Variation::SampleGenotypeFeature;
use Bio::EnsEMBL::Variation::Variation;
use Bio::EnsEMBL::Slice;

our $verbose = 0;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');

my $vdb = $multi->get_DBAdaptor('variation');
my $cdb = $multi->get_DBAdaptor('core');

my $sgfa = $vdb->get_SampleGenotypeFeatureAdaptor();

# test constructor
my $ind = Bio::EnsEMBL::Variation::Individual->new(
  -name => 'test individual',
  -description => 'This is a test individual',
  -gender => 'Male');

my $sample = Bio::EnsEMBL::Variation::Sample->new(
  -name => 'test sample',
  -description => 'This is a test sample',
  -individual => $ind,
);

my $var_id = 1748253;
my $var = Bio::EnsEMBL::Variation::Variation->new(
  -dbID => $var_id,
  -name => 'rs2299222',
  -source => 'dbSNP');


##need a slice
my $pos = 86442404;
my $sa = $cdb->get_SliceAdaptor();
my $slice = $sa->fetch_by_region('chromosome', '7', $pos, $pos, 1);


my $allele1  = 'C';
my $allele2  = 'A';
my @genotype = ($allele1,$allele2);

my $sample_gt_feature = Bio::EnsEMBL::Variation::SampleGenotypeFeature->new(
  -genotype => \@genotype,
  -variation => $var,
  -_variation_id => $var_id,
  -sample => $sample,
  -slice => $slice
);

ok($sample_gt_feature->allele1() eq $allele1, "new - allele 1");
ok($sample_gt_feature->allele2() eq $allele2, "new - allele 2");
ok($sample_gt_feature->genotype_string() eq join('|',@genotype), "new - genotype_string");

ok($sample_gt_feature->sample()->name() eq $sample->name(), "new - sample name");
ok($sample_gt_feature->variation()->name()  eq $var->name(), "new - var name");

my $sample_gt_feature2 = Bio::EnsEMBL::Variation::SampleGenotypeFeature->new_fast({
	'start' => $pos,
	'end' => $pos,
	'strand' => 1,
	'slice' => $slice,
	'genotype' => \@genotype,
	'adaptor' => $sgfa,
	'sample' => $sample,
	'_variation_id' => $var_id
});
			

ok($sample_gt_feature2->allele1() eq $allele1, "new fast - allele 1");
ok($sample_gt_feature2->allele2() eq $allele2, "new fast - allele 2");
ok($sample_gt_feature2->genotype_string() eq join('|',@genotype), "new fast - genotype_string");

ok($sample_gt_feature2->sample()->name() eq $sample->name(), "new fast - sample name");
ok($sample_gt_feature2->variation()->name()  eq $var->name(), "new fast - var name");

done_testing();
