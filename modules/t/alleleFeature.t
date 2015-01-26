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
use FindBin qw($Bin);

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::Individual;

use_ok('Bio::EnsEMBL::Variation::AlleleFeature');
use_ok('Bio::EnsEMBL::Variation::DBSQL::AlleleFeatureAdaptor');

## examples to be added to test-genome-DBs files

=head

my $reg = 'Bio::EnsEMBL::Registry';
$reg->load_all("$Bin/test.ensembl.registry.72");

my $cdba = $reg->get_DBAdaptor('human', 'core');
my $vdba = $reg->get_DBAdaptor('human', 'variation');

my $sa = $cdba->get_SliceAdaptor();
my $slice = $sa->fetch_by_region('chromosome', '1', 1, 200_000);

my $strain_name = 'VENTER';
my $strain_slice = $slice->get_by_strain($strain_name);
my $afs = $strain_slice->get_all_AlleleFeatures_Slice();
is(scalar @$afs, 37, 'Number of AlleleFeatures');

my $hash;
foreach my $af ( @{$afs} ) {
    my $start = $af->start;
    my $end   = $af->end;
    my $strand = $af->strand;
    my $slice = $af->slice;
    my $allele_string = $af->allele_string;
    my $variation_name = $af->variation_name;
    my $variation = $af->variation;
    my $source = $af->source;
    my $individual = $af->individual->name;
    my $consequence_type = join(', ', @{$af->consequence_type});
    my $tvs = $af->get_all_TranscriptVariations;
    my $vf = $af->variation_feature;
    $hash->{$start . '-' . $end . '-' . $strand}->{'allele_string'} = $allele_string;
    $hash->{$start . '-' . $end . '-' . $strand}->{'variation_name'} = $variation_name;
    $hash->{$start . '-' . $end . '-' . $strand}->{'individual'} = $individual;
    $hash->{$start . '-' . $end . '-' . $strand}->{'consequence_type'} = $consequence_type;
}

is($hash->{'29436-29436-1'}->{'allele_string'}, 'G|G', 'Test allele_string');
is($hash->{'29436-29436-1'}->{'consequence_type'}, 'SARA', 'Test consequence_type');
is($hash->{'126113-126113-1'}->{'variation_name'}, 'rs79114531', 'Test variation_name');
is($hash->{'126113-126113-1'}->{'consequence_type'}, 'nc_transcript_variant, intron_variant, downstream_gene_variant, upstream_gene_variant', 'Test consequence_type');
=cut
done_testing();
