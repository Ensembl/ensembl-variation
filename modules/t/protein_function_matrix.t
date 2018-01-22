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

BEGIN {
    use_ok(
        'Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix', 
    );
}

my %preds = (
    polyphen => [
        ['probably damaging', '1'],
        ['probably damaging', '0.931'],
        ['possibly damaging', '0.82'],
        ['possibly damaging', '0.6'],
        ['benign', '0.132'],
        ['benign', '0.12'],
        ['benign', '0'],
        ['unknown', '0'],
    ],
    sift => [
        ['tolerated', '1'],
        ['tolerated', '0.92'],
        ['tolerated', '0.835'],
        ['deleterious', '0.041'],
        ['deleterious', '0.001'],
        ['deleterious', '0'],
    ],
);

# test encoding and decoding prediction as a short

for my $tool (keys %preds) {

    my $m = Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix->new(-analysis => $tool);

    for my $test (@{ $preds{$tool} }) {
        my ($pred, $score) = @$test;
        my $short = $m->prediction_to_short($pred, $score);
        my ($pred2, $score2) =  $m->prediction_from_short($short);
        is($pred, $pred2, "Predictions match: '$pred' vs. '$pred2'");
        is($score, $score2, "Scores match: $score vs. $score2");
    }
}

# test serializing and deserializing a full matrix

my $test_prot = "MACLWY";

my $i = 0;

my $pm = Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix->new(
    -analysis        => 'polyphen',
    -peptide_length  => length($test_prot),
);

my $sm = Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix->new(
    -analysis        => 'sift',
    -peptide_length  => length($test_prot),
);

for my $ref (split //, $test_prot) {
    
    for my $alt (@Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix::ALL_AAS) {
        next if $alt eq $ref;
        $pm->add_prediction($i+1, $alt, @{ $preds{polyphen}->[$i] });
        $sm->add_prediction($i+1, $alt, @{ $preds{sift}->[$i] });
    }

    $i++;
}

my $pph_serialized  = $pm->serialize;
my $sift_serialized = $sm->serialize;

my $pm2 = Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix->new(
    -analysis   => 'polyphen',
    -matrix     => $pph_serialized
);

my $sm2 = Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix->new(
    -analysis   => 'sift',
    -matrix     => $sift_serialized,
);

$i = 0;

for my $pos (1 .. length($test_prot)) {

    for my $aa (@Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix::ALL_AAS) {

        next if $aa eq substr($test_prot, $pos - 1, 1);

        my ($pph_pred, $pph_score)   = $pm2->get_prediction($pos, $aa);
        my ($sift_pred, $sift_score) = $sm2->get_prediction($pos, $aa);
            
        my ($pph_pred2, $pph_score2)   = @{ $preds{polyphen}->[$i] };
        my ($sift_pred2, $sift_score2) = @{ $preds{sift}->[$i] };

        unless( $sift_pred eq $sift_pred2 && $sift_score == $sift_score2 && 
                $pph_pred eq $pph_pred2 && $pph_score == $pph_score2 ) {
            fail("Mismatching results") or die;
        }
    }

    $i++;
}

pass('check serializing and deserializing predictions');

done_testing();

