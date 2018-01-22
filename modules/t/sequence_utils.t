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
use Test::Exception;
use Bio::EnsEMBL::Test::MultiTestDB;

BEGIN {
    use_ok('Bio::EnsEMBL::Variation::Utils::Sequence', qw(sequence_with_ambiguity align_seqs trim_sequences get_matched_variant_alleles));
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



## get_matched_variant_alleles
##############################

# test data integrity checks
throws_ok {get_matched_variant_alleles()} qr/undef/, 'get_matched_variant_alleles - missing a';
throws_ok {get_matched_variant_alleles({})} qr/undef/, 'get_matched_variant_alleles - missing b';

throws_ok {get_matched_variant_alleles([])} qr/expected.+HASH/, 'get_matched_variant_alleles - wrong ref type a';
throws_ok {get_matched_variant_alleles({}, [])} qr/expected.+HASH/, 'get_matched_variant_alleles - wrong ref type b';

throws_ok {get_matched_variant_alleles({}, {})} qr/Missing ref key.+first/, 'get_matched_variant_alleles - missing ref field a';
throws_ok {get_matched_variant_alleles({ref => 'A'}, {})} qr/Missing ref key.+second/, 'get_matched_variant_alleles - missing ref field b';

throws_ok {get_matched_variant_alleles({ref => 'A'}, {ref => 'A'})} qr/Missing alt.+first/, 'get_matched_variant_alleles - missing alts field a';
throws_ok {get_matched_variant_alleles({ref => 'A', alts => ['B']}, {ref => 'A'})} qr/Missing alt.+second/, 'get_matched_variant_alleles - missing alts field b';

throws_ok {get_matched_variant_alleles({ref => 'A', alts => ['B']}, {ref => 'A', alts => ['B']})} qr/Missing pos key.+first/, 'get_matched_variant_alleles - missing pos field a';
throws_ok {get_matched_variant_alleles({ref => 'A', alts => ['B'], pos => 1}, {ref => 'A', alts => ['B']})} qr/Missing pos key.+second/, 'get_matched_variant_alleles - missing pos field b';

# define tests
# {
#   a => input variant A,
#   b => input variant B,
#   r => expected result,
#   d => test description
# }
my @tests = (

  # basic
  {
    a => {ref => 'A', alt => 'G', pos => 1},
    b => {ref => 'A', alt => 'G', pos => 1},
    r => [{
      a_allele => 'G',
      a_index  => 0,
      b_allele => 'G',
      b_index  => 0
    }],
    d => 'basic SNP'
  },
  {
    a => {ref => 'A', alt => '-', pos => 1},
    b => {ref => 'A', alt => '-', pos => 1},
    r => [{
      a_allele => '-',
      a_index  => 0,
      b_allele => '-',
      b_index  => 0
    }],
    d => 'basic del'
  },
  {
    a => {ref => '-', alt => 'A', pos => 1},
    b => {ref => '-', alt => 'A', pos => 1},
    r => [{
      a_allele => 'A',
      a_index  => 0,
      b_allele => 'A',
      b_index  => 0
    }],
    d => 'basic ins'
  },

  # allele_string
  {
    a => {ref => 'A', alt => 'G', pos => 1},
    b => {allele_string => 'A/G', pos => 1},
    r => [{
      a_allele => 'G',
      a_index  => 0,
      b_allele => 'G',
      b_index  => 0
    }],
    d => 'basic SNP - allele_string'
  },
  {
    a => {allele_string => 'A/-', pos => 1},
    b => {ref => 'A', alt => '-', pos => 1},
    r => [{
      a_allele => '-',
      a_index  => 0,
      b_allele => '-',
      b_index  => 0
    }],
    d => 'basic del - allele_string'
  },
  {
    a => {ref => '-', alt => 'A', pos => 1},
    b => {allele_string => '-/A', pos => 1},
    r => [{
      a_allele => 'A',
      a_index  => 0,
      b_allele => 'A',
      b_index  => 0
    }],
    d => 'basic ins - allele_string'
  },


  # multi alleleic
  {
    a => {ref => 'A', alt => 'G', pos => 1},
    b => {ref => 'A', alts => ['C', 'G'], pos => 1},
    r => [{
      a_allele => 'G',
      a_index  => 0,
      b_allele => 'G',
      b_index  => 1
    }],
    d => 'multialleleic SNP 1'
  },
  {
    a => {ref => 'A', alts => ['T', 'G'], pos => 1},
    b => {ref => 'A', alts => ['C', 'G'], pos => 1},
    r => [{
      a_allele => 'G',
      a_index  => 1,
      b_allele => 'G',
      b_index  => 1
    }],
    d => 'multialleleic SNP 2'
  },
  {
    a => {ref => 'A', alts => ['T', 'G'], pos => 1},
    b => {ref => 'A', alts => ['G'], pos => 1},
    r => [{
      a_allele => 'G',
      a_index  => 1,
      b_allele => 'G',
      b_index  => 0
    }],
    d => 'multialleleic SNP 3'
  },

  # rev strand
  {
    a => {ref => 'A', alts => ['G'], pos => 1, strand => -1},
    b => {ref => 'T', alts => ['C'], pos => 1},
    r => [{
      a_allele => 'G',
      a_index  => 0,
      b_allele => 'C',
      b_index  => 0
    }],
    d => 'rev SNP 1'
  },
  {
    a => {ref => 'A', alts => ['G'], pos => 1},
    b => {ref => 'T', alts => ['C'], pos => 1, strand => -1},
    r => [{
      a_allele => 'G',
      a_index  => 0,
      b_allele => 'C',
      b_index  => 0
    }],
    d => 'rev SNP 2'
  },

  # insertions vs VCF-type entries
  {
    a => {ref => '-', alts => ['T'], pos => 2},
    b => {ref => 'A', alts => ['AT'], pos => 1},
    r => [{
      a_allele => 'T',
      a_index  => 0,
      b_allele => 'AT',
      b_index  => 0
    }],
    d => 'ins vs VCF'
  },
  {
    a => {ref => '-', alts => ['TCG'], pos => 2},
    b => {ref => 'A', alts => ['ATCG'], pos => 1},
    r => [{
      a_allele => 'TCG',
      a_index  => 0,
      b_allele => 'ATCG',
      b_index  => 0
    }],
    d => 'ins vs VCF longer'
  },
  {
    a => {ref => '-', alts => ['T', 'TT'], pos => 2},
    b => {ref => 'A', alts => ['AT', 'ATT'], pos => 1},
    r => [
      {
        a_allele => 'T',
        a_index  => 0,
        b_allele => 'AT',
        b_index  => 0
      },
      {
        a_allele => 'TT',
        a_index  => 1,
        b_allele => 'ATT',
        b_index  => 1
      }
    ],
    d => 'ins vs VCF multi 1'
  },
  {
    a => {ref => '-', alts => ['T', 'TT'], pos => 2},
    b => {ref => 'A', alts => ['ATT'], pos => 1},
    r => [{
      a_allele => 'TT',
      a_index  => 1,
      b_allele => 'ATT',
      b_index  => 0
    }],
    d => 'ins vs VCF multi 2'
  },
  {
    a => {ref => '-', alts => ['TT'], pos => 2},
    b => {ref => 'A', alts => ['AT', 'ATT'], pos => 1},
    r => [{
      a_allele => 'TT',
      a_index  => 0,
      b_allele => 'ATT',
      b_index  => 1
    }],
    d => 'ins vs VCF multi 3'
  },
  {
    a => {ref => 'A', alts => ['AT'], pos => 1},
    b => {ref => '-', alts => ['T'], pos => 2},
    r => [{
      a_allele => 'AT',
      a_index  => 0,
      b_allele => 'T',
      b_index  => 0
    }],
    d => 'ins vs VCF other way'
  },

  # deletions vs VCF-type entries
  {
    a => {ref => 'T', alts => ['-'], pos => 2},
    b => {ref => 'AT', alts => ['A'], pos => 1},
    r => [{
      a_allele => '-',
      a_index  => 0,
      b_allele => 'A',
      b_index  => 0
    }],
    d => 'del vs VCF'
  },
  {
    a => {ref => 'TCGT', alts => ['-'], pos => 2},
    b => {ref => 'ATCGT', alts => ['A'], pos => 1},
    r => [{
      a_allele => '-',
      a_index  => 0,
      b_allele => 'A',
      b_index  => 0
    }],
    d => 'del vs VCF longer'
  },
  {
    a => {ref => 'T', alts => ['-'], pos => 2},
    b => {ref => 'AT', alts => ['ATC', 'A'], pos => 1},
    r => [{
      a_allele => '-',
      a_index  => 0,
      b_allele => 'A',
      b_index  => 1
    }],
    d => 'del vs VCF multi 1'
  },

  # trimming
  {
    a => {ref => 'C', alts => ['T'], pos => 3},
    b => {ref => 'ATC', alts => ['ATT'], pos => 1},
    r => [{
      a_allele => 'T',
      a_index  => 0,
      b_allele => 'ATT',
      b_index  => 0
    }],
    d => 'trim SNP 1'
  },
  {
    a => {ref => 'C', alts => ['T'], pos => 3},
    b => {ref => 'CGG', alts => ['TGG'], pos => 3},
    r => [{
      a_allele => 'T',
      a_index  => 0,
      b_allele => 'TGG',
      b_index  => 0
    }],
    d => 'trim SNP 2'
  },
  {
    a => {ref => 'C', alts => ['T'], pos => 3},
    b => {ref => 'ATCGG', alts => ['ATTGG'], pos => 1},
    r => [{
      a_allele => 'T',
      a_index  => 0,
      b_allele => 'ATTGG',
      b_index  => 0
    }],
    d => 'trim SNP 3'
  },
  {
    a => {ref => 'ATC', alts => ['ATT', 'GTC'], pos => 1},
    b => {ref => 'CGG', alts => ['TGG', 'GGG'], pos => 3},
    r => [{
      a_allele => 'ATT',
      a_index  => 0,
      b_allele => 'TGG',
      b_index  => 0
    }],
    d => 'trim SNP 4'
  },

  # horrible mixed things
  {
    a => {ref => 'A', alts => ['-', 'AA'], pos => 2},
    b => {ref => 'TA', alts => ['TAA', 'T'], pos => 1},
    r => [
      {
        a_allele => 'AA',
        a_index  => 1,
        b_allele => 'TAA',
        b_index  => 0
      },{
        a_allele => '-',
        a_index  => 0,
        b_allele => 'T',
        b_index  => 1
      },
    ],
    d => 'ins and del 1'
  },
  {
    a => {ref => 'CA', alts => ['-', 'CACA'], pos => 2},
    b => {ref => 'GCA', alts => ['GCACA', 'G'], pos => 1},
    r => [
      {
        a_allele => 'CACA',
        a_index  => 1,
        b_allele => 'GCACA',
        b_index  => 0
      },{
        a_allele => '-',
        a_index  => 0,
        b_allele => 'G',
        b_index  => 1
      },
    ],
    d => 'ins and del 2'
  },
  {
    a => {ref => 'A', alts => ['-', 'G', 'AA'], pos => 2},
    b => {ref => 'TA', alts => ['TG', 'TAA', 'T'], pos => 1},
    r => [
      {
        a_allele => 'G',
        a_index  => 1,
        b_allele => 'TG',
        b_index  => 0
      },
      {
        a_allele => 'AA',
        a_index  => 2,
        b_allele => 'TAA',
        b_index  => 1
      },
      {
        a_allele => '-',
        a_index  => 0,
        b_allele => 'T',
        b_index  => 2
      },
    ],
    d => 'ins, del and SNP'
  },
  {
    a => {ref => 'T', alts => ['-', 'C', 'AA', 'TT'], pos => 2, strand => -1},
    b => {ref => 'TA', alts => ['TG', 'TAA', 'T'], pos => 1},
    r => [
      {
        a_allele => 'C',
        a_index  => 1,
        b_allele => 'TG',
        b_index  => 0
      },
      {
        a_allele => 'TT',
        a_index  => 3,
        b_allele => 'TAA',
        b_index  => 1
      },
      {
        a_allele => '-',
        a_index  => 0,
        b_allele => 'T',
        b_index  => 2
      },
    ],
    d => 'ins, del and SNP (reverse strand)'
  },
  {
    a => {ref => 'A', alts => ['-', 'G', 'AA'], pos => 2},
    b => {ref => 'TA', alts => ['TG', 'CAA', 'T'], pos => 1},
    r => [
      {
        a_allele => 'G',
        a_index  => 1,
        b_allele => 'TG',
        b_index  => 0
      },
      {
        a_allele => '-',
        a_index  => 0,
        b_allele => 'T',
        b_index  => 2
      },
    ],
    d => 'ins (mismatched 1st base), del and SNP'
  },

  # some real world examples from ExAC
  {
    a => {allele_string => '-/GGCGGC', pos => 3257703},
    b => {allele_string => 'AGGC/A/AGGCGGCGGC', pos => 3257699},
    r => [{
      a_allele => 'GGCGGC',
      a_index  => 0,
      b_allele => 'AGGCGGCGGC',
      b_index  => 1
    }],
    d => '-/GGCGGC vs AGGC/A/AGGCGGCGGC'
  },
  {
    a => {allele_string => 'A/-', pos => 112253194},
    b => {allele_string => 'CAAA/CAAAA/CA/CAA/CAAAAA/C/CAAAAAA', pos => 112253191},
    r => [{
      a_allele => '-',
      a_index  => 0,
      b_allele => 'CAA',
      b_index  => 2
    }],
    d => 'A/- vs CAAA/CAAAA/CA/CAA/CAAAAA/C/CAAAAAA'
  },
  {
    a => {allele_string => '-/ACAC/ACACAC', pos => 33485790},
    b => {allele_string => 'AACACACACAC/AAC/AACACACAC/AACACACACACACACAC/A/AACACACACACACAC/AACACAC', pos => 33485779},
    r => [
      {
        a_allele => 'ACACAC',
        a_index  => 1,
        b_allele => 'AACACACACACACACAC',
        b_index  => 2
      },
      {
        a_allele => 'ACAC',
        a_index  => 0,
        b_allele => 'AACACACACACACAC',
        b_index  => 4
      },
    ],
    d => '-/ACAC/ACACAC vs AACACACACAC/AAC/AACACACAC/AACACACACACACACAC/A/AACACACACACACAC/AACACAC'
  },
  {
    a => {allele_string => '-/TG/TGTG/TGTGTG', pos => 231700228},
    b => {allele_string => 'TTGTGTG/TTGTGTGTGTG/TTGTG/TTGTGTGTG/TTG/TTGTGTGTGTGTG/T', pos => 231700221},
    r => [
      {
        a_allele => 'TGTG',
        a_index  => 1,
        b_allele => 'TTGTGTGTGTG',
        b_index  => 0
      },
      {
        a_allele => 'TG',
        a_index  => 0,
        b_allele => 'TTGTGTGTG',
        b_index  => 2
      },
      {
        a_allele => 'TGTGTG',
        a_index  => 2,
        b_allele => 'TTGTGTGTGTGTG',
        b_index  => 4
      },
    ],
    d => '-/TG/TGTG/TGTGTG vs TTGTGTG/TTGTGTGTGTG/TTGTG/TTGTGTGTG/TTG/TTGTGTGTGTGTG/T'
  },
  {
    a => {allele_string => 'A/AA/-', pos => 179117285},
    b => {allele_string => 'TAA/TA/TAAA/TAAAA/T', pos => 179117284},
    r => [
      {
        'a_index' => 1,
        'a_allele' => '-',
        'b_allele' => 'TA',
        'b_index' => 0
      },
      {
        'a_index' => 0,
        'a_allele' => 'AA',
        'b_allele' => 'TAAA',
        'b_index' => 1
      }
    ],
    d => 'A/AA/- vs TAA/TA/TAAA/TAAAA/T'
  },
);

# run tests
foreach my $test(@tests) {
  is_deeply(
    get_matched_variant_alleles($test->{a}, $test->{b}),
    $test->{r},
    'get_matched_variant_alleles - '.($test->{d} || 'misc')
  );
}


## sequence with ambiguity
my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');
my $vdb = $multi->get_DBAdaptor('variation');
my $cdb = $multi->get_DBAdaptor('core');

my $seq = sequence_with_ambiguity($cdb, $vdb, 9, 22124503, 22126503);
my $exp = 'NWNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNMNNNNNNNNNNNNNNNNNNNNNNNNRNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNSNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNRNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNMNNNNNNNNNNNNRNNNNNNNRNNNNNNNNNNNNNNNRNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNYNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNSNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNWNNNNMNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNYNNNNNNNNNNNNNNNNNNNNNNNNNNSNNNNNNNNNNNNNMNNNNNNNNNNNNNNNRNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNRNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNKNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNYNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNYNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNRNNNNNNNNNNYNNNSNNNNNNNNNNNNNNNNYNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNMNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNRNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN';

ok($seq->isa('Bio::EnsEMBL::Slice'), "sequence_with_ambiguity isa slice");
is($seq->seq, $exp, "sequence_with_ambiguity seq");



done_testing();

