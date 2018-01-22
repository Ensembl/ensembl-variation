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
    use_ok('Bio::EnsEMBL::Variation::Utils::Sequence', qw(variation_class));
}

my %tests = (
    'A/T'                   => 'snp',
    'G/C/A'                 => 'snp',
    'A'                     => 'het',
    'cnv'                   => 'cnv',
    'CNV_PROBE'             => 'cnv probe',
    'HGMD_MUTATION'         => 'hgmd_mutation',
    'AA/TTT'                => 'substitution',
    '-/A'                   => 'in-del',
    '-/(LARGEINSERTION)'    => 'named',
    '-/(2345 BP INSERTION)' => 'named',
    '-/INSERTION'           => 'named',
    '-/(LARGEDELETION)'     => 'named', # dbSNP stylee
    '(1657 BP DELETION)/-'  => 'named', # COSMIC stylee
    'ATTAGC/-'              => 'in-del',
    'A/-'                   => 'in-del',
    'DELETION/-'            => 'named',
    '(LARGEDELETION)/-/AT'  => 'mixed',
    'AGCG/-/A'              => 'mixed',
    'A/-/T'                 => 'mixed',
    'A/Y/-'                 => 'mixed',
    '(CA)1/-/(CA)12'        => 'microsat',
    '-/TGTG/(TG)10/TG(11)'  => 'mixed',
);

for my $allele_string (keys %tests) {
    my $expected = $tests{$allele_string};
    is(variation_class($allele_string, 0), $expected, "$allele_string => $expected") ;
}

is(variation_class('A/G', 1), 'somatic_snv', "somatic");

done_testing();

