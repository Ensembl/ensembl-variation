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

use Bio::SimpleAlign;
use Bio::LocatableSeq;

BEGIN {
    use_ok('Bio::EnsEMBL::Variation::Utils::ComparaUtils', qw(dump_alignment_for_polyphen dump_alignment_for_sift));
}

# we should turn this:

my %input = (
    QUERY   => '-ABC--DEFG-H---', 
    S1      => '---IJKLMN--OPQR',
    S2      => 'ST------UVWHYZ-',
);

# into this:

my %output = (
    QUERY   => 'ABCDEFGH',
    S1      => '--ILMN-O',
    S2      => 'T----UVH',
);

# create a test alignment

my $sa = Bio::SimpleAlign->new;

for my $k (sort keys %input) {

    $sa->add_seq(
        Bio::LocatableSeq->new(
            -SEQ    => $input{$k},
            -START  => 1,
            -END    => length($input{$k}),
            -ID     => $k,
            -STRAND => 0
        )
    );
}

# and ungap it (polyphen excludes the query seq, while sift includes it)

my $sift_align = Bio::EnsEMBL::Variation::Utils::ComparaUtils::_ungap_alignment($sa, 'QUERY', 1);
my ($query_seq, $pph_align) = Bio::EnsEMBL::Variation::Utils::ComparaUtils::_ungap_alignment($sa, 'QUERY', 0);

# we just check that the sift alignment looks good

for my $n ($sift_align->each_seq) {
    is($n->seq, $output{$n->display_id}, $input{$n->display_id}.' => '.$output{$n->display_id});
}

# check the percent_id;

is(
    Bio::EnsEMBL::Variation::Utils::ComparaUtils::_percent_id($query_seq, $pph_align->get_seq_by_pos(1)), 
    0, 
    "S1 percent_id is correct"
);

is(
    Bio::EnsEMBL::Variation::Utils::ComparaUtils::_percent_id($query_seq, $pph_align->get_seq_by_pos(2)), 
    0.125,
    "S2 percent_id is correct"
);

# check we create alignment files in the expected format

my $sift_expected = <<SIFT;
>QUERY/1-8
ABCDEFGH
>S1/1-5
--ILMN-O
>S2/1-4
T----UVH
SIFT

my $pph_expected = <<PPH;
CLUSTAL QUERY (2)

S2                                                                    T----UVH
S1                                                                    --ILMN-O
PPH

{
    no warnings qw(redefine);
    
    # redefine the _get_ungapped_alignment subroutine to return our test alignment 

    sub Bio::EnsEMBL::Variation::Utils::ComparaUtils::_get_ungapped_alignment {
        return wantarray ? ($query_seq, $pph_align) : $sift_align;
    }
}

undef $/;

my $pph_file = 'pph.aln';
my $sift_file = 'sift.fa';

dump_alignment_for_polyphen('QUERY', $pph_file);
ok(-e $pph_file, "polyphen file created");
open my $PPH, "<$pph_file";
is(<$PPH>, $pph_expected, "polyphen file contents look correct");
`rm $pph_file`;

dump_alignment_for_sift('QUERY', $sift_file);
ok(-e $sift_file, "sift file created");
open my $SIFT, "<$sift_file";
is(<$SIFT>, $sift_expected, "sift file contents look correct");
`rm $sift_file`;

done_testing();

