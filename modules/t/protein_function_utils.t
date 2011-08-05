use strict;
use warnings;

use Test::More;

BEGIN {
    use_ok(
        'Bio::EnsEMBL::Variation::Utils::ProteinFunctionUtils', 
        qw(
            @ALL_AAS
            $AA_LOOKUP
            $NO_PREDICTION
            $HEADER
            prediction_from_string 
            prediction_to_short
            prediction_from_short
            compress_prediction_string
            expand_prediction_string
            compute_offset
        )
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
    for my $test (@{ $preds{$tool} }) {
        my ($pred, $score) = @$test;
        my $short = prediction_to_short($tool, $pred, $score);
        my ($pred2, $score2) =  prediction_from_short($tool, $short);
        is($pred, $pred2, "Predictions match: '$pred' vs. '$pred2'");
        is($score, $score2, "Scores match: $score vs. $score2");
    }
}

my $sift_string;
my $pph_string;

my $test_prot = "MACLWY";

my $i = 0;

for my $ref (split //, $test_prot) {
    for my $alt (@ALL_AAS) {
        if ($alt eq $ref) {
            $sift_string .= $NO_PREDICTION;
            $pph_string  .= $NO_PREDICTION;
        }
        else {
            $sift_string .= prediction_to_short('sift', @{ $preds{sift}->[$i] });
            $pph_string  .= prediction_to_short('polyphen', @{ $preds{polyphen}->[$i] });
        }
   }

    $i++;
}

my $sift_string_zipped = compress_prediction_string($sift_string);
my $pph_string_zipped  = compress_prediction_string($pph_string);

ok(length($sift_string_zipped) < length($sift_string), "sift string compresses");
ok(length($pph_string_zipped) < length($pph_string), "polyphen string compresses");

my $sift_string_unzipped = expand_prediction_string($sift_string_zipped);
my $pph_string_unzipped  = expand_prediction_string($pph_string_zipped);

is($sift_string_unzipped, $HEADER.$sift_string, "sift string decompresses ok");
is($pph_string_unzipped, $HEADER.$pph_string, "polyphen string decompresses ok") or die;

$i = 0;

for my $pos (1 .. length($test_prot)) {

    for my $aa (@ALL_AAS) {

        next if $aa eq substr($test_prot, $pos - 1, 1);

        my ($sift_pred, $sift_score) = prediction_from_string('sift', $sift_string_unzipped, $pos, $aa);
        my ($pph_pred, $pph_score)   = prediction_from_string('polyphen', $pph_string_unzipped, $pos, $aa);

        my ($sift_pred2, $sift_score2) = @{ $preds{sift}->[$i] };
        my ($pph_pred2, $pph_score2) = @{ $preds{polyphen}->[$i] };

        unless($sift_pred eq $sift_pred2 && $sift_score == $sift_score2 && $pph_pred eq $pph_pred2 && $pph_score == $pph_score2) {
            fail("Mismatching results") or die;
        }
    }

    $i++;
}

pass('check reading predictions from a string');

done_testing();

