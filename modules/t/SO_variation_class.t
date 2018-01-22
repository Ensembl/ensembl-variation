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
    use_ok('Bio::EnsEMBL::Variation::Utils::Sequence', qw(SO_variation_class));
}

my %tests = (
    'A/T'                   => 'SNV',
    'A/W'                   => 'SNV',
    'M/K'                   => 'SNV',
    'G/C/A'                 => 'SNV',
    'G/C/W'                 => 'SNV',
    'AA/TT'                 => 'substitution',
    'A/YY'                  => 'indel',
    'GA/TTT/ACCCC'          => 'indel',
    '-/A'                   => 'insertion',
    '-/TAAG'                => 'insertion',
    '-/(LARGEINSERTION)'    => 'insertion',
    '-/(2345 BP INSERTION)' => 'insertion',
    '-/W'                   => 'insertion',
    '-/AA/ATGCG'            => 'insertion',
    '(508 BP INSERTION)'    => 'insertion',
    '-/INSERTION'           => 'insertion',
    '-/(LARGEDELETION)'     => 'deletion', # dbSNP stylee
    '(1657 BP DELETION)/-'  => 'deletion', # COSMIC stylee
    'ATTAGC/-'              => 'deletion',
    'A/-'                   => 'deletion',
    'YY/-'                  => 'deletion',
    '(LARGEDELETION)'       => 'deletion',
    'DELETION/-'            => 'deletion',
    '(LARGEDELETION)/-/AT'  => 'sequence_alteration',
    'AGCG/-/A'              => 'sequence_alteration',
    'A/-/T'                 => 'sequence_alteration',
    'A/Y/-'                 => 'sequence_alteration',
    '(CA)1/-/(CA)12'        => 'tandem_repeat',
    '(CAG)8/(CAG)9'         => 'tandem_repeat',
    '-/TGTG/(TG)10/TG(11)'  => 'tandem_repeat',
    '-/(RY)7/(RY)8'         => 'tandem_repeat',
);

for my $allele_string (keys %tests) {
    my $expected = $tests{$allele_string};
    is(SO_variation_class($allele_string, 1), $expected, "$allele_string => $expected") ;
    if ($expected =~ /insertion|deletion/) {
        is(SO_variation_class($allele_string, 0), 'indel', "$allele_string => indel if ref_correct = 0") ;
    }
}

done_testing();

