# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2019] EMBL-European Bioinformatics Institute
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
    use_ok('Bio::EnsEMBL::Variation::Utils::Sequence', qw(SO_variation_class));
}

my %tests = (
    'A/T'                   => 'SNV',
    'A/W'                   => 'SNV',
    'M/K'                   => 'SNV',
    'G/C/A'                 => 'SNV',
    'G/C/W'                 => 'SNV',
    'G/*/C'                 => 'SNV',
    'AA/TT'                 => 'substitution',
    'A/YY'                  => 'indel',
    'GA/TTT/ACCCC'          => 'indel',
    'CTCT/CTCTCT'           => 'insertion',
    'CTCT/CT'               => 'deletion',
    'CTCT/CTCTCT/CT'        => 'indel',
    'CT/-'                  => 'deletion',
    'CT/CTCT'               => 'insertion',
    'CT/-/CTCT'             => 'indel', # Not a sequence alteration
    'C/-/CC'                => 'indel',
    'TCT/T/TCTG'            => 'indel',
    'TCT/-/TCTG'            => 'indel',
    '-/A'                   => 'insertion',
    '-/TAAG'                => 'insertion',
    '-/(LARGEINSERTION)'    => 'insertion',
    '-/(2345 BP INSERTION)' => 'insertion',
    '-/W'                   => 'insertion',
    '-/AA/ATGCG'            => 'insertion',
    '(508 BP INSERTION)'    => 'insertion',
    '-/INSERTION'           => 'insertion',
    'T/TT'                  => 'insertion',
    'TT/TTTT'               => 'insertion',
    'TT/TTTT/TTTTT'         => 'insertion',
    '-/(LARGEDELETION)'     => 'deletion', # dbSNP style
    '(1657 BP DELETION)/-'  => 'deletion', # COSMIC style
    'ATTAGC/-'              => 'deletion',
    'A/-'                   => 'deletion',
    'YY/-'                  => 'deletion',
    '(LARGEDELETION)'       => 'deletion',
    'DELETION/-'            => 'deletion',
    'TT/T'                  => 'deletion',
    'TTT/TT/T'              => 'deletion',
    'TT/-/T'                => 'deletion',
    'CTT/C'                 => 'deletion',
    'TA/T'                  => 'deletion',
    'TCT/CT'                => 'deletion',
    'AGCG/-/A'              => 'deletion',
    '(LARGEDELETION)/-/AT'  => 'sequence_alteration',
    'A/-/T'                 => 'sequence_alteration',
    'A/Y/-'                 => 'sequence_alteration',
    '(CA)1/-/(CA)12'        => 'tandem_repeat',
    '(CAG)8/(CAG)9'         => 'tandem_repeat',
    '-/TGTG/(TG)10/TG(11)'  => 'tandem_repeat',
    '-/(RY)7/(RY)8'         => 'tandem_repeat',
    'TT/T/TTT'              => 'indel',     # dbSNP v2 style indel
    'TT/TTT/TTTT'           => 'insertion', # dbSNP v2 style indel
    'TTT/T/TT'              => 'deletion',  # dbSNP v2 stype indel
);

for my $allele_string (sort keys %tests) {
    my $expected = $tests{$allele_string};
    is(SO_variation_class($allele_string, 1), $expected, "$allele_string => $expected") ;
    if ($expected =~ /insertion|deletion/) {
        is(SO_variation_class($allele_string, 0), 'indel', "$allele_string => indel if ref_correct = 0") ;
    }
}

done_testing();
