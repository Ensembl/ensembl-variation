# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2017] EMBL-European Bioinformatics Institute
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
use Test::Exception;
use Bio::EnsEMBL::Test::MultiTestDB;

BEGIN {
    use_ok('Bio::EnsEMBL::Variation::Utils::Sequence', qw(sequence_with_ambiguity align_seqs trim_sequences));
}



## align seqs
my $s1 = 'ACGTACGT';
my $s2 = 'AGGTACGT';

my $aligned = align_seqs($s1, $s2);
is($s1, $aligned->[0], "align_seqs - sub 1");
is($s2, $aligned->[1], "align_seqs - sub 2");

$s2 = 'ATACGT';
$aligned = align_seqs($s1, $s2);
is($s1, $aligned->[0], "align_seqs - del 1");
is($aligned->[1], 'A--TACGT', "align_seqs - del 2");

$s2 = 'AGGTACG';
$aligned = align_seqs($s1, $s2);
is($s1, $aligned->[0], "align_seqs - end del 1");
is($aligned->[1], 'AGGTACG-', "align_seqs - end del 2");

$s2 = 'CGTACGT';
$aligned = align_seqs($s1, $s2);
is($s1, $aligned->[0], "align_seqs - start del 1");
is($aligned->[1], '-CGTACGT', "align_seqs - start del 2");



## trim_sequences
#################

is_deeply(
  trim_sequences(qw(A B)),
  [qw(A B 0 0 0)],
  'trim_sequences - no change'
);

is_deeply(
  trim_sequences(qw(CA CB)),
  [qw(A B 1 1 1)],
  'trim_sequences - beginning'
);

is_deeply(
  trim_sequences(qw(AC BC)),
  [qw(A B 0 0 1)],
  'trim_sequences - end'
);

is_deeply(
  trim_sequences(qw(DAC DBC)),
  [qw(A B 1 1 1)],
  'trim_sequences - both'
);

is_deeply(
  trim_sequences(qw(FOOABAR FOOBBAR)),
  [qw(A B 3 3 1)],
  'trim_sequences - both long'
);

is_deeply(
  trim_sequences(qw(DAC DBC 10)),
  [qw(A B 11 11 1)],
  'trim_sequences - coords'
);

is_deeply(
  trim_sequences(qw(ATTT AT 10 13 0 0)),
  ['TT', '', 12, 13, 1],
  'trim_sequences - trim from left first'
);

is_deeply(
  trim_sequences(qw(ATTT AT 10 13 0 1)),
  ['TT', '', 11, 12, 1],
  'trim_sequences - trim from right first'
);

is_deeply(
  trim_sequences(qw(ATTT AT 10 13 1)),
  ['TT', '-', 12, 13, 1],
  'trim_sequences - empty_to_dash'
);

throws_ok {trim_sequences(undef, 'A')} qr/Missing reference or alternate sequence/, 'trim_sequences - no ref';
throws_ok {trim_sequences('A')} qr/Missing reference or alternate sequence/, 'trim_sequences - no alt';
throws_ok {trim_sequences()} qr/Missing reference or alternate sequence/, 'trim_sequences - no both';



## sequence with ambiguity
my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');
my $vdb = $multi->get_DBAdaptor('variation');
my $cdb = $multi->get_DBAdaptor('core');

my $seq = sequence_with_ambiguity($cdb, $vdb, 9, 22124503, 22126503);
my $exp = 'NWNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNMNNNNNNNNNNNNNNNNNNNNNNNNRNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNSNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNRNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNMNNNNNNNNNNNNRNNNNNNNRNNNNNNNNNNNNNNNRNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNYNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNSNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNWNNNNMNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNYNNNNNNNNNNNNNNNNNNNNNNNNNNSNNNNNNNNNNNNNMNNNNNNNNNNNNNNNRNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNRNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNKNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNYNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNYNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNRNNNNNNNNNNYNNNSNNNNNNNNNNNNNNNNYNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNMNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNRNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN';

ok($seq->isa('Bio::EnsEMBL::Slice'), "sequence_with_ambiguity isa slice");
is($seq->seq, $exp, "sequence_with_ambiguity seq");



done_testing();

