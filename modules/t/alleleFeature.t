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
use FindBin qw($Bin);

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Variation::Sample;

use_ok('Bio::EnsEMBL::Variation::AlleleFeature');
use_ok('Bio::EnsEMBL::Variation::DBSQL::AlleleFeatureAdaptor');

## examples to be added to test-genome-DBs files

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');
my $vdb = $multi->get_DBAdaptor('variation');
my $cdb = $multi->get_DBAdaptor('core');

my $chr = '9';
my $sa = $cdb->get_SliceAdaptor();
my $slice = $sa->fetch_by_region('chromosome',$chr,22124503,22126503);

my $strain_slice_adaptor = $vdb->get_StrainSliceAdaptor;

my $strain_name = '1000GENOMES:phase_1:NA06984';

my $strain_slice = $strain_slice_adaptor->get_by_strain_Slice($strain_name, $slice);
my $seq = $strain_slice->seq();
is(substr($seq, 1, 1), "W", "apply_edit (via sequence)");

my $afs = $strain_slice->get_all_AlleleFeatures_Slice();
is(scalar @$afs, 27, 'Number of AlleleFeatures');

my $v_name = 'rs1333047';
my $cons   = 'downstream_gene_variant';


my $af1 = $afs->[0];
print "AF: ".$af1->ref_allele_string."\n";
ok($af1->start == 2,                       'allele - start');
ok($af1->end == 2,                         'allele - end');
ok($af1->strand == 1,                      'allele - strand');
ok($af1->slice->seq_region_name eq $chr,   'allele - slice');
ok($af1->allele_string eq 'A|T',           'allele - allele_string');
ok($af1->variation_name eq $v_name,        'allele - variation_name');
ok($af1->variation->name eq $v_name,       'allele - variation');
ok($af1->source eq 'dbSNP',                'allele - source');
ok($af1->sample->name eq $strain_name,     'allele - sample');
ok($af1->consequence_type->[0] eq $cons,   'allele - consequence_type');
ok($af1->display_consequence eq $cons,     'allele - display_consequence');
ok($af1->ref_allele_string eq 'N',         'allele - ref_allele_string');
ok($af1->length == 1,                      'allele - length');
ok($af1->length_diff == 0,                 'allele - length_diff');

# test get all OverlapConsequences
my $overlap_cons = $af1->get_all_OverlapConsequences();
ok($overlap_cons->[0]->SO_term eq $cons,   'allele - get_all_OverlapConsequences');

# test most severe OverlapConsequence
my $msc = $af1->most_severe_OverlapConsequence();
ok($msc->SO_term eq $cons, 'allele - most_severe_OverlapConsequence');

# test get_all_TranscriptVariations (no data available in test db)
my $tvs = $af1->get_all_TranscriptVariations;
ok($tvs->[0]->transcript_stable_id eq 'ENST00000422420', 'allele - get_all_TranscriptVariations');

# test variation feature
my $vf = $af1->variation_feature;
ok($vf->variation_name eq $v_name, 'allele - variation_feature');

# test get all sources
my $sources = $af1->get_all_sources();
ok($sources->[0] eq 'dbSNP', 'allele - get_all_sources');


my $hash;
foreach my $af (@$afs) {
    my $start = $af->start;
    my $end   = $af->end;
    my $strand = $af->strand;
    my $slice = $af->slice;
    my $allele_string = $af->allele_string;
    my $variation_name = $af->variation_name;
    my $variation = $af->variation;
    my $source = $af->source;
    my $sample = $af->sample->name;
    my $consequence_type = join(', ', @{$af->consequence_type});
    my $vf = $af->variation_feature;
    $hash->{$start . '-' . $end . '-' . $strand}->{'allele_string'} = $allele_string;
    $hash->{$start . '-' . $end . '-' . $strand}->{'variation_name'} = $variation_name;
    $hash->{$start . '-' . $end . '-' . $strand}->{'sample'} = $sample;
    $hash->{$start . '-' . $end . '-' . $strand}->{'consequence_type'} = $consequence_type;
}

print "\n";
is($hash->{'1411-1411-1'}->{'allele_string'}, 'C|T', 'Test allele_string');
is($hash->{'1001-1001-1'}->{'consequence_type'}, 'downstream_gene_variant', 'Test consequence_type');
is($hash->{'1083-1083-1'}->{'variation_name'}, 'rs73650063', 'Test variation_name');


done_testing();
