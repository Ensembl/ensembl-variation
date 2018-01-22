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
use Data::Dumper;
use FindBin qw($Bin);
use File::Path qw(remove_tree);;

use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Variation::Utils::FastaSequence qw(setup_fasta clear_fasta_cache revert_fasta);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::CoordSystem;
use Bio::EnsEMBL::Slice;

my $gz_fasta = "$Bin\/testdata/vep-cache/homo_sapiens/78_GRCh38/test.fa.gz";
my $fasta = $gz_fasta;
$fasta =~ s/\.gz//;
`gzip -dc $gz_fasta > $fasta` if(-e "$gz_fasta");

my ($db, $slice1, $slice2, $seq1, $seq2);

$db = setup_fasta(-FASTA => $fasta);
ok($db, "basic");
ok($db->isa('Bio::DB::HTS::Faidx') || $db->isa('Bio::DB::Fasta'), "isa");

# different way of calling
$db = setup_fasta($fasta);
ok($db, "single arg");

# throws
throws_ok(sub {setup_fasta()}, qr/No FASTA file specified/, 'throws - No FASTA file specified');
throws_ok(sub {setup_fasta(-FASTA => '/does/not/exist')}, qr/not found/, 'throws - not found');
throws_ok(sub {setup_fasta(-FASTA => $fasta, -TYPE => 'invalid')}, qr/Unrecognised index type/, 'throws - Unrecognised index type');
# throws_ok(sub {setup_fasta(-FASTA => $gz_fasta, -TYPE => 'Bio::DB::Fasta')}, qr/Cannot index bgzipped FASTA file with Bio::DB::Fasta/, 'throws - Cannot index bgzipped FASTA file with Bio::DB::Fasta');


my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');
my $cdb = $multi->get_DBAdaptor('core');
my $sa = $cdb->get_SliceAdaptor();

$db = setup_fasta(-FASTA => $fasta, -SYNONYMS => {foo => {21 => 1}});

$slice1 = $sa->fetch_by_region('chromosome', 21, 25606454, 25606454);
is($slice1->seq, 'G', "single bp slice");

is($slice1->expand(0, 5)->seq, 'GCACAA', "expand single bp slice 3'");
is($slice1->expand(5, 0)->seq, 'GCATGG', "expand single bp slice 5'");
is($slice1->expand(5, 5)->seq, 'GCATGGCACAA', "expand single bp slice 5' and 3'");

$seq1 = $slice1->expand(5, 5)->seq;
reverse_comp(\$seq1);

clear_fasta_cache();
is($slice1->expand(5, 5)->invert->seq, $seq1, "expand single bp slice 5' and 3' and invert after cache clear");


# do same tests again, starting with reverse strand slice
clear_fasta_cache();

$slice1 = $sa->fetch_by_region('chromosome', 21, 25606454, 25606454, -1);
is($slice1->seq, 'C', "rev - single bp slice");

is($slice1->expand(0, 5)->seq, 'CCATGC', "rev - expand single bp slice 3'");
is($slice1->expand(5, 0)->seq, 'TTGTGC', "rev - expand single bp slice 5'");
is($slice1->expand(5, 5)->seq, 'TTGTGCCATGC', "rev - expand single bp slice 5' and 3'");

$seq1 = $slice1->expand(5, 5)->seq;
reverse_comp(\$seq1);

clear_fasta_cache();
is($slice1->expand(5, 5)->invert->seq, $seq1, "rev - expand single bp slice 5' and 3' and invert after cache clear");


# test subseq
clear_fasta_cache();
$slice1 = $sa->fetch_by_region('chromosome', 21, 25606450, 25606460);
is($slice1->subseq, 'CATGGCACAAC', 'subseq - no params');
is($slice1->subseq(1, 5), 'CATGG', 'subseq 1');
is($slice1->subseq(3, 5), 'TGG', 'subseq 2');
is($slice1->subseq(3, 5, -1), 'CCA', 'subseq rev');
is($slice1->subseq(5, 4), '', 'subseq e > s');
is($slice1->subseq(-1, 5), 'CGCATGG', "subseq overlap 5'");


# test going off ends
$slice1 = $sa->fetch_by_region('chromosome', 21, 1, 10);
is($slice1->seq, 'N' x 10, "start of chrom");
is($slice1->expand(10, 0)->seq, 'N' x 20, "expand beyond start of chrom");


# test synonyms
clear_fasta_cache();
$slice1 = Bio::EnsEMBL::Slice->new(
  -COORD_SYSTEM      => Bio::EnsEMBL::CoordSystem->new(-NAME => 'chromosome', -RANK => 1),
  -START             => 25606454,
  -END               => 25606454,
  -SEQ_REGION_NAME   => 'foo',
  -SEQ_REGION_LENGTH => $slice1->length
);
is($slice1->seq, 'G', "synonym");

$slice2 = Bio::EnsEMBL::Slice->new(
  -COORD_SYSTEM      => Bio::EnsEMBL::CoordSystem->new(-NAME => 'chromosome', -RANK => 1),
  -START             => 25606454,
  -END               => 25606454,
  -SEQ_REGION_NAME   => 'chrfoo',
  -SEQ_REGION_LENGTH => $slice1->length
);
is($slice2->seq, 'G', "implied chr synonym");

# remember to revert!!!
revert_fasta();

# and clean up
unlink($fasta);
unlink("$fasta\.index");
unlink("$fasta\.fai");

done_testing();
