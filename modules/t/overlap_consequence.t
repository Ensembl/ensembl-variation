# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2023] EMBL-European Bioinformatics Institute
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
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Test::MultiTestDB;
my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');
my $cdba = $multi->get_DBAdaptor('core');
my $vdba = $multi->get_DBAdaptor('variation');

my $vfa = $vdba->get_VariationFeatureAdaptor;
my $slice_adaptor = $cdba->get_SliceAdaptor;
my $chrom = 13;
my $slice = $slice_adaptor->fetch_by_region('chromosome', $chrom);
my $vf = Bio::EnsEMBL::Variation::VariationFeature->new_fast({
    start          => 51483899,
    end            => 51484001,
    allele_string  =>  "GCCCTGCGGCGCGCGATGAGCACCTACCACTAGCGGAGCCGCGAGGGAGAGGCCGCGGCCCCTTCCCGTTGCCTGCGGCCACCGGCCGGCATTCAGAGCCCCT/-",
    strand         => 1,
    map_weight     => 1,
    slice          => $slice,
    adaptor        => $vfa,
});
my $consequence_hash = {};
my $vfos = $vf->get_all_VariationFeatureOverlaps;
foreach my $vfo (@{$vfos}) {
  my $vfoas = $vfo->get_all_alternate_VariationFeatureOverlapAlleles;
  foreach my $vfoa (@$vfoas) {
    my $consequences = $vfoa->get_all_OverlapConsequences;
    foreach my $consequence (@$consequences) {
      $consequence_hash->{$consequence->SO_term} = 1;
    }
  }
}
my $expected = '5_prime_UTR_variant,non_coding_transcript_exon_variant,regulatory_region_ablation,regulatory_region_variant';
my $got = join(',', sort keys %$consequence_hash);
ok($got eq $expected, "All overlap consequences for variation feature overlaps");

$vf = Bio::EnsEMBL::Variation::VariationFeature->new_fast({
    start          => 51483899,
    end            => 51484001,
    allele_string  =>  "GCCCTGCGGCGCGCGATGAGCACCTACCACTAGCGGAGCCGCGAGGGAGAGGCCGCGGCCCCTTCCCGTTGCCTGCGGCCACCGGCCGGCATTCAGAGCCCCT/-",
    strand         => 1,
    map_weight     => 1,
    slice          => $slice,
    adaptor        => $vfa,
});
$consequence_hash = {};
my $rfvs = $vf->get_all_RegulatoryFeatureVariations;
foreach my $rfv (@$rfvs) {
  my $rfvas = $rfv->get_all_alternate_RegulatoryFeatureVariationAlleles;
  foreach my $rfva (@$rfvas) {
    my $consequences = $rfva->get_all_OverlapConsequences;
    foreach my $consequence (@$consequences) {
      $consequence_hash->{$consequence->SO_term} = 1;
    }
  }
}

$expected = 'regulatory_region_ablation,regulatory_region_variant';
$got = join(',', sort keys %$consequence_hash);
ok($got eq $expected, "All overlap consequences for regulatory variation features");

done_testing();
