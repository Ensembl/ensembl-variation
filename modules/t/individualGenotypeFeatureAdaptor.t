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

our $verbose = 0;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');

my $vdba = $multi->get_DBAdaptor('variation');
my $cdba = $multi->get_DBAdaptor('core');

my $va   = $vdba->get_VariationAdaptor();
my $igfa = $vdba->get_IndividualGenotypeFeatureAdaptor();
my $inda = $vdba->get_IndividualAdaptor();
my $popa = $vdba->get_PopulationAdaptor();
my $sa   = $cdba->get_SliceAdaptor();

ok($igfa && $igfa->isa('Bio::EnsEMBL::Variation::DBSQL::IndividualGenotypeFeatureAdaptor'), 'isa IndividualGenotypeFeatureAdaptor');


my $ind_name = '1000GENOMES:phase_1:NA07048';
my $ind = $inda->fetch_all_by_name($ind_name)->[0];

my $pop_name = '1000GENOMES:phase_1_CEU';
my $pop = $popa->fetch_by_name($pop_name);

my $var_name = 'rs1333047';#'rs7569578';
my $variation = $va->fetch_by_name($var_name);

# fetch_all_by_Variation with Individual
my $igfs = $igfa->fetch_all_by_Variation($variation, $ind);
ok($igfs->[0]->individual->name eq $ind_name, 'fetch_all_by_Variation - with individual');

# fetch_all_by_Slice with Individual
my $slice = $sa->fetch_by_region('chromosome','9',22124504,22124504);
my $igfs2 = $igfa->fetch_all_by_Slice($slice, $ind);
ok($igfs2->[0]->variation->name eq $var_name, 'fetch_all_by_Slice - with individual');

# fetch_all_by_Slice with Population
my $igfs3 = $igfa->fetch_all_by_Slice($slice, $pop);
my @igf_ind = grep {$_->individual->name eq $ind_name} @$igfs3;
ok($igf_ind[0]->individual->name eq $ind_name, 'fetch_all_by_Slice - with population');

done_testing();
