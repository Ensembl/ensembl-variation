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
use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Variation::Individual;

use_ok('Bio::EnsEMBL::Variation::AlleleFeature');
use_ok('Bio::EnsEMBL::Variation::DBSQL::AlleleFeatureAdaptor');

## examples to be added to test-genome-DBs files

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');
my $vdb = $multi->get_DBAdaptor('variation');
my $cdb = $multi->get_DBAdaptor('core');

my $sa = $cdb->get_SliceAdaptor();
my $slice = $sa->fetch_by_region('chromosome','9',22124503,22126503);

my $strain_name = '1000GENOMES:phase_1:NA06984';

my $strain_slice = $slice->get_by_strain($strain_name);
my $seq = $strain_slice->seq();
is(substr($seq, 1, 1), "W", "apply_edit (via sequence)");

my $afs = $strain_slice->get_all_AlleleFeatures_Slice();
is(scalar @$afs, 27, 'Number of AlleleFeatures');

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

is($hash->{'1411-1411-1'}->{'allele_string'}, 'C|T', 'Test allele_string');
is($hash->{'1001-1001-1'}->{'consequence_type'}, 'downstream_gene_variant', 'Test consequence_type');
is($hash->{'1083-1083-1'}->{'variation_name'}, 'rs73650063', 'Test variation_name');

my $af = $afs->[0];
is($af->display_consequence, 'downstream_gene_variant', "display_consequence");

done_testing();
