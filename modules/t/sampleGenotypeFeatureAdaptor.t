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

our $verbose = 0;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');

my $vdba = $multi->get_DBAdaptor('variation');
my $cdba = $multi->get_DBAdaptor('core');

my $va   = $vdba->get_VariationAdaptor();
my $sgfa = $vdba->get_SampleGenotypeFeatureAdaptor();
my $sample_adpt = $vdba->get_SampleAdaptor();
my $pa = $vdba->get_PopulationAdaptor();
my $sa   = $cdba->get_SliceAdaptor();

ok($sgfa && $sgfa->isa('Bio::EnsEMBL::Variation::DBSQL::SampleGenotypeFeatureAdaptor'), 'isa SampleGenotypeFeatureAdaptor');

my $sample_name = '1000GENOMES:phase_1:NA07048';
my $sample = $sample_adpt->fetch_all_by_name($sample_name)->[0];

my $pop_name = '1000GENOMES:phase_1_CEU';
my $pop = $pa->fetch_by_name($pop_name);

my $var_name = 'rs1333047';#'rs7569578';
my $variation = $va->fetch_by_name($var_name);

# fetch_all_by_Variation with Sample
my $sgfs = $sgfa->fetch_all_by_Variation($variation, $sample);
ok($sgfs->[0]->sample->name eq $sample_name, 'fetch_all_by_Variation - with sample');

throws_ok { $sgfa->fetch_all_by_Variation('variation'); } qr/Bio::EnsEMBL::Variation::Variation argument expected/, 'Throw on wrong argument for fetch_all_by_Variation';
warns_like {
  $sgfa->fetch_all_by_Variation(Bio::EnsEMBL::Variation::Variation->new(-name => 'variation_name'));
} qr/Cannot retrieve genotypes for variation without set dbID/, 'Warn if wrong argument for fetch_all_by_Variation';

# fetch_all_by_Slice with Sample
my $slice = $sa->fetch_by_region('chromosome','9',22124504,22124504);
my $sgfs2 = $sgfa->fetch_all_by_Slice($slice, $sample);
ok($sgfs2->[0]->variation->name eq $var_name, 'fetch_all_by_Slice - with sample');

# fetch_all_by_Slice with Population
my $sgfs3 = $sgfa->fetch_all_by_Slice($slice, $pop);
my @sgf_sample = grep {$_->sample->name eq $sample_name} @$sgfs3;
ok($sgf_sample[0]->sample->name eq $sample_name, 'fetch_all_by_Slice - with population');


done_testing();
