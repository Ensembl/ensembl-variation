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
my $mfva = $vdb->get_MotifFeatureVariationAdaptor;

my $sa = $cdb->get_SliceAdaptor();
my $rfa = $fdb->get_RegulatoryFeatureAdaptor();
my $mfa = $fdb->get_MotifFeatureAdaptor();

my $slice = $sa->fetch_by_region('chromosome', '7');
my $v = Bio::EnsEMBL::Variation::Variation->new(
  -name => 'rs140471675_mf_test',
  -_source_id => 1,
  -is_somatic => 0,
);
$va->store($v);
$v = $va->fetch_by_name('rs140471675_mf_test');
my $vf = Bio::EnsEMBL::Variation::VariationFeature->new(
  -start   => 151409224,
  -end     => 151409224,
  -strand  => 1,
  -slice   => $slice,
  -allele_string => 'G/A',
  -variation_name => 'rs140471675_mf_test',
  -map_weight  => 1,
  -_source_id => 1,
  -is_somatic => 0,
  -variation => $v,
);
$vfa->store($vf);

my $feature_id = 939053;
my $motif_feature = $mfa->fetch_by_dbID($feature_id) or die "Failed to fetch MotifFeature for id: $feature_id";
my $rf = $rfa->fetch_all_by_attribute_feature($motif_feature)->[0];
my $mfv = Bio::EnsEMBL::Variation::MotifFeatureVariation->new(
    -motif_feature      => $motif_feature,
    -variation_feature  => $vf,
    -adaptor            => $mfva,
    -disambiguate_single_nucleotide_alleles => 0,
);
if ($mfv && (scalar(@{$mfv->consequence_type}) > 0) ) {
    $mfva->store($mfv, $rf);
}

my $consequence_types = join(',', sort @{$mfv->consequence_type});
ok($consequence_types eq 'TF_binding_site_variant', 'Print consequence types for motif_feature_variation');

my $motif_name = $mfv->motif_name;
ok($motif_name eq 'Max:MA0058.1', 'Compare motif name');
my $feature_stable_id = $mfv->feature_stable_id;
ok($feature_stable_id eq 'ENSR00000636355', 'Compare feature stable id');

my $alt_MotifFeatureVariationAlleles = $mfv->get_all_alternate_MotifFeatureVariationAlleles;

$feature_id = 769450;
$motif_feature = $mfa->fetch_by_dbID($feature_id) or die "Failed to fetch MotifFeature for id: $feature_id";

my $mfvs = $mfva->fetch_all_by_MotifFeatures([$motif_feature]);
ok(scalar @$mfvs == 1, 'fetch_all_by_MotifFeatures');

my $mfvs1 = $mfva->fetch_all_somatic_by_MotifFeatures([$motif_feature]);
ok(scalar @$mfvs1 == 0, 'fetch_all_somatic_by_MotifFeatures');

my $mfvs2 = $mfva->fetch_all_by_MotifFeatures_SO_terms([$motif_feature], ['sequence_variant']);
ok(scalar @$mfvs2 == 0, 'fetch_all_by_MotifFeatures_SO_terms');

my $mfvs3 = $mfva->fetch_all_somatic_by_MotifFeatures_SO_terms([$motif_feature], ['sequence_variant']);
ok(scalar @$mfvs3 == 0, 'fetch_all_somatic_by_MotifFeatures_SO_terms');

my $var4 = $va->fetch_by_name('rs372423729');
my $vf4 = $vfa->fetch_all_by_Variation($var4)->[0];

#my $mfvs4 = $mfva->fetch_all_by_VariationFeatures_SO_terms([$vf4], [$motif_feature], ['sequence_variant']);
#ok(scalar @$mfvs4 == 0, 'fetch_all_by_VariationFeatures_SO_terms');

#my $count = $mfva->count_all_by_VariationFeatures_SO_terms([$vf4], [$motif_feature], ['sequence_variant']);
#ok($count == 0, 'count_all_by_VariationFeatures_SO_terms');

done_testing();
