# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2025] EMBL-European Bioinformatics Institute
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

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');

my $vdb = $multi->get_DBAdaptor('variation');
my $cdb = $multi->get_DBAdaptor('core');

my $svf_adaptor   = $vdb->get_StructuralVariationFeatureAdaptor;
my $slice_adaptor = $cdb->get_SliceAdaptor;

use_ok('Bio::EnsEMBL::Variation::StructuralVariation');
use_ok('Bio::EnsEMBL::Variation::StructuralVariationFeature');
use_ok('Bio::EnsEMBL::Variation::StructuralVariationOverlap');
use_ok('Bio::EnsEMBL::Variation::StructuralVariationOverlapAllele');

# test copy number variant
my $dbID = 4509635;
my $outer_start = 7803891;
my $inner_start = 7805991;
my $inner_end = 7823440;
my $outer_end = 7825340;
my $var_name = 'esv93078';
my $chr = '8';
my $is_somatic = 0;
my $sv_length = $outer_end-$outer_start+1;
my $SO_term = 'copy_number_variant';

my $svf = Bio::EnsEMBL::Variation::StructuralVariationFeature->new
  (-adaptor     => $svf_adaptor,
   -outer_start => $outer_start,
   -start       => $outer_start,
   -inner_start => $inner_start,
   -inner_end   => $inner_end,
   -end         => $outer_end,
   -outer_end   => $outer_end,
   -strand      => 1,
   -slice       => $slice_adaptor->fetch_by_region('chromosome', $chr),
   -variation_name => $var_name,
   -is_somatic => $is_somatic,
   -length => $sv_length,
   -class_SO_term => $SO_term,
);

sub _create_feature {
  my $obj = shift;
  my $chr = $obj->{chr} || $obj->seq_region_name;
  return Bio::EnsEMBL::Feature->new(
    -seqname => $chr,
    -start   => $obj->{start} - 100,
    -end     => $obj->{end} + 100,
    -strand  => 0,
    -slice   => $slice_adaptor->fetch_by_region('chromosome', $chr),
  );
}

my $feature = _create_feature($svf);
my $svo = Bio::EnsEMBL::Variation::StructuralVariationOverlap->new(
  -feature                      => $feature,
  -structural_variation_feature => $svf,
  -no_transfer                  => 1,
);
is($svo->{breakends}, undef, "cnv -> no breakends");
my $svoas = $svo->get_all_StructuralVariationOverlapAlleles();
is($svoas->[0]->allele_number, 1, "cnv -> allele_number 1");
is($svoas->[0]->structural_variation_overlap, $svo, "cnv -> svo");

# test chromosomal breakpoint
$svf = Bio::EnsEMBL::Variation::StructuralVariationFeature->new(
  -adaptor         => $svf_adaptor,
  -outer_start     => 7803891,
  -start           => 7803891,
  -inner_start     => 7803891,
  -inner_end       => 7803891,
  -end             => 7803891,
  -outer_end       => 7803891,
  -seq_region_name => '8',
  -slice           => $slice_adaptor->fetch_by_region('chromosome', 8),
  -strand          => 1,
  -allele_string   => "A[4:66578[/T[7:3433[",
  -class_SO_term   => 'chromosome_breakpoint',
);
$svf->{chr} = $svf->seq_region_name;

# create SVOverlap object
my $feat0 = _create_feature($svf);
$svo = Bio::EnsEMBL::Variation::StructuralVariationOverlap->new(
  -feature                      => $feat0,
  -structural_variation_feature => $svf,
  -no_transfer                  => 1,
);
$svoas = $svo->get_all_StructuralVariationOverlapAlleles();
is($svoas->[0]->symbolic_allele, 'N.', "svoa -> ref symbolic allele");

# check if breakends were parsed after calling SVOverlap
is_deeply($svf->{breakends}->[0]->{chr},       4, "svf -> breakend 1 chr");
is_deeply($svf->{breakends}->[0]->{start}, 66578, "svf -> breakend 1 start");
is_deeply($svf->{breakends}->[1]->{chr},       7, "svf -> breakend 2 chr");
is_deeply($svf->{breakends}->[1]->{start},  3433, "svf -> breakend 2 start");

# test if variation/breakends are close to specific features
sub _close_to_feature {
  Bio::EnsEMBL::Variation::StructuralVariationOverlap::_close_to_feature(@_);
}

my $bnd1  = $svf->{breakends}->[0];
my $bnd2  = $svf->{breakends}->[1];
is(_close_to_feature($svf,  $feat0), 1, "svf -> near feature 0");
is(_close_to_feature($bnd1, $feat0), 0, "svf breakend 1 -> far from feature 0");
is(_close_to_feature($bnd2, $feat0), 0, "svf breakend 2 -> far from feature 0");

my $feat1 = _create_feature($bnd1);
is(_close_to_feature($svf,  $feat1), 0, "svf -> far from feature 1");
is(_close_to_feature($bnd1, $feat1), 1, "svf breakend 1 -> near feature 1");
is(_close_to_feature($bnd2, $feat1), 0, "svf breakend 2 -> far from feature 1");

my $feat2 = _create_feature($bnd2);
is(_close_to_feature($svf,  $feat2), 0, "svf -> far from feature 2");
is(_close_to_feature($bnd1, $feat2), 0, "svf breakend 1 -> far from feature 2");
is(_close_to_feature($bnd2, $feat2), 1, "svf breakend 2 -> near feature 2");

# create SVOverlap object with feat1
$svo = Bio::EnsEMBL::Variation::StructuralVariationOverlap->new(
  -feature                      => $feat1,
  -structural_variation_feature => $svf,
  -no_transfer                  => 1,
);
$svoas = $svo->get_all_StructuralVariationOverlapAlleles();
is($svoas->[0]->symbolic_allele, 'A[4:66578[', "svoa -> breakend 1 symbolic allele");


# test single breakend
$svf = Bio::EnsEMBL::Variation::StructuralVariationFeature->new(
  -adaptor         => $svf_adaptor,
  -outer_start     => 7803891,
  -start           => 7803891,
  -inner_start     => 7803891,
  -inner_end       => 7803891,
  -end             => 7803891,
  -outer_end       => 7803891,
  -seq_region_name => '8',
  -slice           => $slice_adaptor->fetch_by_region('chromosome', 8),
  -strand          => 1,
  -allele_string   => "TCG.",
  -class_SO_term   => 'chromosome_breakpoint',
);
$svf->{chr} = $svf->seq_region_name;

# create SVOverlap object
$feat0 = _create_feature($svf);
$svo = Bio::EnsEMBL::Variation::StructuralVariationOverlap->new(
  -feature                      => $feat0,
  -structural_variation_feature => $svf,
  -no_transfer                  => 1,
);
$svoas = $svo->get_all_StructuralVariationOverlapAlleles();
is($svoas->[0]->base_variation_feature->{string}, 'TCG.', "svoa -> single breakend");

# test if variation/breakends are close to specific features
is(_close_to_feature($svf,  $feat0), 1, "svf -> near feature 0");

done_testing();
