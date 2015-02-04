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

my $va = $vdba->get_VariationAdaptor();
my $igta  = $vdba->get_IndividualGenotypeAdaptor();

ok($va && $va->isa('Bio::EnsEMBL::Variation::DBSQL::VariationAdaptor'), 'isa variation_adaptor');
ok($igta && $igta->isa('Bio::EnsEMBL::Variation::DBSQL::IndividualGenotypeAdaptor'), 'isa individual_genotype_adaptor');

# fetch_all_by_Variation
my $variation_name = 'rs7569578';
my $variation = $va->fetch_by_name($variation_name);

my $igts = $igta->fetch_all_by_Variation($variation);
ok(scalar @$igts == 3, 'number of IGTs for variation');

my @filtered_igts = grep {$_->individual->name eq 'NA18635' } @$igts;
ok(scalar @filtered_igts == 1, 'number of IGTs for individual');

my $igt = $filtered_igts[0];
ok($igt->allele1 eq 'A', 'allele1');
ok($igt->allele2 eq 'T', 'allele2');


# store
$igt = $igts->[0];
my $i = $igt->individual->adaptor->fetch_all_by_name('GM18507')->[0];
$igt->individual($i);

$igta->store([$igt], 1);

delete($igta->{_cache});
$igts = $igta->fetch_all_by_Variation($variation);
($igt) = grep {$_->individual->name eq 'GM18507'} @$igts;

ok($igt && $igt->genotype_string eq 'T|T', "store merged");

$igt = $igts->[0];

ok($igta->store([$igt]), "store unmerged");
ok($igta->store_uncompressed([$igt]), "store uncompressed");

# fetch_all_by_Slice
#my $sa = $cdba->get_SliceAdaptor;
#my $slice = $sa->fetch_by_region('chromosome', 13, 54950000, 55000000);

done_testing();
