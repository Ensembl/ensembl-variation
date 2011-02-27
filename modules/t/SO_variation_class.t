use strict;
use warnings;

use Test::More;

BEGIN {
    use_ok('Bio::EnsEMBL::Variation::Utils::Sequence', qw(SO_variation_class));
}

my %tests = (
    'A/T'                   => 'SNV',
    'G/C/A'                 => 'SNV',
    'AA/TTT'                => 'substitution',
    'GA/TTT/ACCCC'          => 'substitution',
    '-/A'                   => 'insertion',
    '-/TAAG'                => 'insertion',
    '-/(LARGEINSERTION)'    => 'insertion',
    '-/(2345 BP INSERTION)' => 'insertion',
    '-/AA/ATGCG'            => 'insertion',
    '(508 BP INSERTION)'    => 'insertion',
    '-/INSERTION'           => 'insertion',
    '-/(LARGEDELETION)'     => 'deletion', # dbSNP stylee
    '(1657 BP DELETION)/-'  => 'deletion', # COSMIC stylee
    'ATTAGC/-'              => 'deletion',
    'A/-'                   => 'deletion',
    '(LARGEDELETION)'       => 'deletion',
    'DELETION/-'            => 'deletion',
    '(LARGEDELETION)/-/AT'  => 'sequence_alteration',
    'AGCG/-/A'              => 'sequence_alteration',
    'A/-/T'                 => 'sequence_alteration',
    '(CA)1/-/(CA)12'        => 'tandem_repeat',
    '(CAG)8/(CAG)9'         => 'tandem_repeat',
    '-/TGTG/(TG)10/TG(11)'  => 'tandem_repeat',
);

for my $allele_string (keys %tests) {
    my $expected = $tests{$allele_string};
    is(SO_variation_class($allele_string, 1), $expected, "$allele_string => $expected") ;
    if ($expected =~ /insertion|deletion/) {
        is(SO_variation_class($allele_string, 0), 'indel', "$allele_string => indel if ref_correct = 0") ;
    }
}

done_testing();

