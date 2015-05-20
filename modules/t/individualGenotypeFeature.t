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
use Bio::EnsEMBL::Variation::IndividualGenotypeFeature;
use Bio::EnsEMBL::Variation::Variation;
use Bio::EnsEMBL::Slice;

our $verbose = 0;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');

my $vdb = $multi->get_DBAdaptor('variation');
my $cdb = $multi->get_DBAdaptor('core');


my $igfa = $vdb->get_IndividualGenotypeFeatureAdaptor();

# test constructor
my $ind = Bio::EnsEMBL::Variation::Individual->new(
  -name => 'test individual',
  -description => 'This is a test individual',
  -gender => 'Male');

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

my $ind_gt_feature = Bio::EnsEMBL::Variation::IndividualGenotypeFeature->new(
  -genotype      => \@genotype,
  -variation     => $var,
  -_variation_id => $var_id,
  -individual    => $ind,
  -slice         => $slice
);

ok($ind_gt_feature->allele1() eq $allele1, "new - allele 1");

done_testing();
