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
my $sgta  = $vdba->get_SampleGenotypeAdaptor();

ok($va && $va->isa('Bio::EnsEMBL::Variation::DBSQL::VariationAdaptor'), 'isa variation_adaptor');
ok($sgta && $sgta->isa('Bio::EnsEMBL::Variation::DBSQL::SampleGenotypeAdaptor'), 'isa sample_genotype_adaptor');

# fetch_all_by_Variation
my $variation_name = 'rs7569578';
my $variation = $va->fetch_by_name($variation_name);

my $sgts = $sgta->fetch_all_by_Variation($variation);
ok(scalar @$sgts == 3, 'number of SGTs for variation');

my @filtered_sgts = grep {$_->sample->name eq 'NA18635' } @$sgts;
ok(scalar @filtered_sgts == 1, 'number of SGTs for sample');

my $sgt = $filtered_sgts[0];
ok($sgt->allele1 eq 'A', 'allele1');
ok($sgt->allele2 eq 'T', 'allele2');


# store
$sgt = $sgts->[0];
my $s = $sgt->sample->adaptor->fetch_all_by_name('GM18507')->[0];
$sgt->sample($s);

$sgta->store([$sgt], 1);

delete($sgta->{_cache});
$sgts = $sgta->fetch_all_by_Variation($variation);
($sgt) = grep {$_->sample->name eq 'GM18507'} @$sgts;

ok($sgt && $sgt->genotype_string eq 'T|T', "store merged");

$sgt = $sgts->[0];

ok($sgta->store([$sgt]), "store unmerged");
ok($sgta->store_uncompressed([$sgt]), "store uncompressed");

# fetch_all_by_Slice
#my $sa = $cdba->get_SliceAdaptor;
#my $slice = $sa->fetch_by_region('chromosome', 13, 54950000, 55000000);

done_testing();
