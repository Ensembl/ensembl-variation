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
use Bio::EnsEMBL::Registry;

our $verbose = 0;
my $omulti = Bio::EnsEMBL::Test::MultiTestDB->new('multi');
my $odb = $omulti->get_DBAdaptor('ontology');
Bio::EnsEMBL::Registry->add_db($omulti, 'ontology', $odb);

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');
my $vdb = $multi->get_DBAdaptor('variation');
my $cdb = $multi->get_DBAdaptor('core');
my $fdb = $multi->get_DBAdaptor('funcgen');
$vdb->dnadb($cdb);

my $va   = $vdb->get_VariationAdaptor();
my $vfa  = $vdb->get_VariationFeatureAdaptor();
my $rfva = $vdb->get_RegulatoryFeatureVariationAdaptor();

my $sa = $cdb->get_SliceAdaptor();

my $rfa = $fdb->get_RegulatoryFeatureAdaptor();

my $rf = $rfa->fetch_by_stable_id('ENSR00000000637');

my $variation_name = 'rs187207343';
my $variation = $va->fetch_by_name($variation_name);
my $vfs = $vfa->fetch_all_by_Variation($variation);
my $vf = $vfs->[0];

my $rfvs = $rfva->fetch_all_by_RegulatoryFeatures([$rf]);
ok(scalar @$rfvs == 1, 'count fetch_all_by_RegulatoryFeatures');

my $rfv = $rfvs->[0];
$rf = $rfv->regulatory_feature;
ok($rf && $rf->isa('Bio::EnsEMBL::Funcgen::RegulatoryFeature'), 'isa regulatory_feature');
my $rf_stable_id = $rfv->regulatory_feature_stable_id;
ok($rf_stable_id eq 'ENSR00000000637', 'regulator_feature stable_id');

my $reference_rfva = $rfv->get_reference_RegulatoryFeatureVariationAllele;
my $ref_allele_string = $reference_rfva->allele_string;
ok($ref_allele_string eq 'C', 'Ref allele string');

my $alt_rfvas = $rfv->get_all_alternate_RegulatoryFeatureVariationAlleles;
ok(scalar @$alt_rfvas == 1, 'count alternate RegulatoryFeatureVariationAlleles');
my $alt_rfva = $alt_rfvas->[0];
my $alt_allele_string = $alt_rfva->allele_string;
ok($alt_allele_string eq 'C/T', 'Alt allele string');

$rfvs = $rfva->fetch_all_somatic_by_RegulatoryFeatures([$rf]);
ok(scalar @$rfvs == 0, 'count fetch_all_somatic_by_RegulatoryFeatures');

#$rfvs = $rfva->fetch_all_by_RegulatoryFeatures_SO_terms([$rf], ['']);

#$rfvs = $rfva->fetch_all_somatic_by_RegulatoryFeatures_SO_terms([$rf]);

$rfvs = $rfva->fetch_all_by_VariationFeatures([$vf]);
ok(scalar @$rfvs == 1, 'count fetch_all_by_VariationFeatures');

#$rfvs = $rfva->fetch_all_by_VariationFeatures_SO_terms([$vf], []);
#my $count = $rfva->count_all_by_VariationFeatures_SO_terms([$vf]);

$rfv = Bio::EnsEMBL::Variation::RegulatoryFeatureVariation->new(
  -regulatory_feature => $rf,
  -variation_feature  => $vf,
  -adaptor            => $rfva,
  -disambiguate_single_nucleotide_alleles => 0,
);
my $consequence_types = join(',', sort @{$rfv->consequence_type});
ok($consequence_types eq 'regulatory_region_variant', 'Print consequence types for regulatory_feature_variation');

my $slice = $sa->fetch_by_region('chromosome', '7');
my $v = $va->fetch_by_name('rs140471675');
if (!$v) {
  $v = Bio::EnsEMBL::Variation::Variation->new(
    -name => 'rs140471675',
    -_source_id => 1,
    -is_somatic => 0,
  );
  $va->store($v);
  $v = $va->fetch_by_name('rs140471675');
}

$vfs = $vfa->fetch_all_by_Variation($v);
$vf = $vfs->[0];
=begin
if (!scalar @$vfs) {
  $vf = Bio::EnsEMBL::Variation::VariationFeature->new(
    -start   => 151409224,
    -end     => 151409224,
    -strand  => 1,
    -slice   => $slice,
    -allele_string => 'G/A',
    -variation_name => 'rs140471675',
    -map_weight  => 1,
    -_source_id => 1,
    -is_somatic => 0,
    -variation => $v,
  );
  $vfa->store($vf);
} else {
  $vf = $vfs->[0];
}
=end
=cut

$rf = $rfa->fetch_by_stable_id('ENSR00000636355');
$rfv = Bio::EnsEMBL::Variation::RegulatoryFeatureVariation->new(
  -regulatory_feature => $rf,
  -variation_feature  => $vf,
  -adaptor            => $rfva,
  -disambiguate_single_nucleotide_alleles => 0,
);

$consequence_types = join(',', sort @{$rfv->consequence_type});
ok($consequence_types eq 'regulatory_region_variant', 'Print consequence types for regulatory_feature_variation');

ok($rfva->store($rfv), 'store');

my $var1 = $va->fetch_by_name('rs187207343');
my $vf1 = $vfa->fetch_all_by_Variation($var1)->[0];
my $rf1 = $rfa->fetch_by_stable_id('ENSR00000000637');
my $rfvs1 = $rfva->fetch_all_by_VariationFeatures_SO_terms([$vf1], [$rf1], ['sequence_variant']);
ok(scalar @$rfvs1 == 0, 'fetch_all_by_VariationFeatures_SO_terms');
my $count = $rfva->count_all_by_VariationFeatures_SO_terms([$vf1], [$rf1], ['sequence_variant']);
ok($count == 0, 'count_all_by_VariationFeatures_SO_terms');
my $rfvs2 = $rfva->fetch_all_by_RegulatoryFeatures_SO_terms([$rf1], ['sequence_variant']);
ok(scalar @$rfvs2 == 0, 'fetch_all_by_RegulatoryFeatures_SO_terms');
my $rfvs3 = $rfva->fetch_all_somatic_by_RegulatoryFeatures_SO_terms([$rf1], ['sequence_variant']);
ok(scalar @$rfvs3 == 0, 'fetch_all_somatic_by_RegulatoryFeatures_SO_terms');
my $rfvs4 = $rfva->fetch_all_by_VariationFeatures([$vf1]);
is(scalar @$rfvs, 1, 'fetch_all_by_VariationFeatures');

done_testing();
