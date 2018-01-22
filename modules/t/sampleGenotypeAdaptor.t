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
use Data::Dumper;

our $verbose = 0;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');

my $vdba = $multi->get_DBAdaptor('variation');
my $cdba = $multi->get_DBAdaptor('core');

my $va = $vdba->get_VariationAdaptor();
my $sgta = $vdba->get_SampleGenotypeAdaptor();
my $sa = $vdba->get_SampleAdaptor();
my $pa = $vdba->get_PopulationAdaptor();

ok($va && $va->isa('Bio::EnsEMBL::Variation::DBSQL::VariationAdaptor'), 'isa variation_adaptor');
ok($sgta && $sgta->isa('Bio::EnsEMBL::Variation::DBSQL::SampleGenotypeAdaptor'), 'isa sample_genotype_adaptor');
ok($sa && $sa->isa('Bio::EnsEMBL::Variation::DBSQL::SampleAdaptor'), 'isa sample_adaptor');
ok($pa && $pa->isa('Bio::EnsEMBL::Variation::DBSQL::PopulationAdaptor'), 'isa population_adaptor');

# fetch_all_by_Variation
my $variation_name = 'rs7569578';
my $variation = $va->fetch_by_name($variation_name);

my $sample_name = 'J. CRAIG VENTER';
my $sample = $sa->fetch_all_by_name($sample_name)->[0];

my $population_name = 'YUSUKE:POP_J';
my $population = $pa->fetch_by_name($population_name);

my $sgts = $sgta->fetch_all_by_Variation($variation, $population);
ok(scalar @$sgts == 0, 'number of SGTs for variation and population');

$sgts = $sgta->fetch_all_by_Variation($variation, $sample);
ok(scalar @$sgts == 1, 'number of SGTs for variation and sample');

throws_ok { $sgta->fetch_all_by_Variation($variation, $variation); } qr/Argument supplied is not of type Bio::EnsEMBL::Variation::Sample or Bio::EnsEMBL::Variation::Population/, 'Throw on wrong argument for fetch_all_by_Variation';

$sgts = $sgta->fetch_all_by_Variation($variation);
ok(scalar @$sgts == 3, 'number of SGTs for variation');

throws_ok { $sgta->fetch_all_by_Variation('variation'); } qr/Variation argument expected/, 'Throw on wrong argument for fetch_all_by_Variation';
warns_like {
  $sgta->fetch_all_by_Variation(Bio::EnsEMBL::Variation::Variation->new(-name => 'variation_name'));
} qr/Cannot retrieve genotypes for variation without dbID/, 'Throw on wrong argument for fetch_all_by_Variation';

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
my $slice_adaptor = $cdba->get_SliceAdaptor;
my $slice = $slice_adaptor->fetch_by_region('chromosome', 13, 54950000, 55000000);
$sgts = $sgta->fetch_all_by_Slice($slice);
ok(scalar @$sgts == 0, 'fetch_all_by_Slice');

# fetch_all_by_Variation_dbID
$sgts = $sgta->fetch_all_by_Variation_dbID($variation->dbID);
ok(scalar @$sgts == 3, 'fetch_all_by_Variation_dbID');

done_testing();
