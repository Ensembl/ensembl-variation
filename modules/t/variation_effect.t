# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2025] EMBL-European Bioinformatics Institute
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

use Data::Dumper;

use FindBin qw($Bin);

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::VariationFeature;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor ;
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::Variation::Utils::VariationEffect;
use Bio::EnsEMBL::Test::MultiTestDB;

BEGIN {
    use_ok('Bio::EnsEMBL::Variation::TranscriptVariation');
}

our $DEBUG = 0;

## test inline C

# test overlap
# coords: f1_start, f1_end, f2_start, f2_end
my @coords = (
  {
    expected => 1,
    coords   => [  1,  1,  1,  1 ],
  },
  {
    expected => 1,
    coords   => [  1,  2,  1,  1 ],
  },
  {
    expected => 1,
    coords   => [  1,  1,  1,  2 ],
  },
  {
    expected => 0,
    coords   => [  1,  1,  2,  2 ],
  },
  {
    expected => 0,
    coords   => [  2,  2,  1,  1 ],
  },
  {
    expected => 1,
    coords   => [  1,  3,  2,  4 ],
  },
  {
    expected => 1,
    coords   => [  2,  4,  1,  3 ],
  },
  {
    expected => 1,
    coords   => [  1,  4,  2,  3 ],
  },
  {
    expected => 1,
    coords   => [  2,  3,  1,  4 ],
  },
  {
    expected => 1,
    coords   => [ -1,  4,  0,  3 ],
  },
);

foreach my $c(@coords) {
  my $p_res   = Bio::EnsEMBL::Variation::Utils::VariationEffect::overlap_perl(@{$c->{coords}}) || 0;
  my $c_res   = Bio::EnsEMBL::Variation::Utils::VariationEffect::overlap(@{$c->{coords}}) || 0;
  my $exp     = $c->{expected};
  my $comment = $c->{comment} || "coords: ".join(", ", @{$c->{coords}});
  
  is($p_res, $exp, "perl overlap $comment");
  is($c_res, $exp, "C    overlap $comment");
}

# test _intron_overlap
# coords: vf_start, vf_end, intron_start, intron_end, bool insertion
@coords = (

  # 3,1 upstream
  {
    expected => 1,
    coords   => [  7,  7, 10, 20, 0 ],
    comment  => "just within 3,1 upstream"
  },
  {
    expected => 1,
    coords   => [  8,  8, 10, 20, 0 ],
    comment  => "compeletely within 3,1 upstream"
  },
  {
    expected => 1,
    coords   => [  9, 10, 10, 20, 0 ],
    comment  => "part 5' overlap 3,1 upstream"
  },
  {
    expected => 1,
    coords   => [  5,  7, 10, 20, 0 ],
    comment  => "part 3' overlap 3,1 upstream"
  },
  {
    expected => 0,
    coords   => [  6,  6, 10, 20, 0 ],
    comment  => "just 5' of 3,1 upstream"
  },
  {
    expected => 0,
    coords   => [ 10, 10, 10, 20, 0 ],
    comment  => "just 3' of 3,1 upstream"
  },
  
  # 1,3 downstream
  {
    expected => 1,
    coords   => [ 33, 33, 10, 30, 0 ],
    comment  => "just within 1,3 downstream"
  },
  {
    expected => 1,
    coords   => [ 32, 32, 10, 30, 0 ],
    comment  => "compeletely within 1,3 downstream"
  },
  {
    expected => 1,
    coords   => [ 30, 31, 10, 30, 0 ],
    comment  => "part 5' overlap 1,3 downstream"
  },
  {
    expected => 1,
    coords   => [ 33, 34, 10, 30, 0 ],
    comment  => "part 3' overlap 1,3 downstream"
  },
  {
    expected => 0,
    coords   => [ 30, 30, 10, 30, 0 ],
    comment  => "just 5' of 1,3 downstream"
  },
  {
    expected => 0,
    coords   => [ 34, 34, 10, 30, 0 ],
    comment  => "just 3' of 1,3 downstream"
  },

  # 2,7 inside 3'
  {
    expected => 1,
    coords   => [ 12, 12, 10, 30, 0 ],
    comment  => "just within 2,7 inside 3'"
  },
  {
    expected => 1,
    coords   => [ 14, 16, 10, 30, 0 ],
    comment  => "compeletely within 2,7 inside 3'"
  },
  {
    expected => 1,
    coords   => [ 12, 14, 10, 30, 0 ],
    comment  => "part 5' overlap 2,7 inside 3'"
  },
  {
    expected => 1,
    coords   => [ 16, 19, 10, 30, 0 ],
    comment  => "part 3' overlap 2,7 inside 3'"
  },
  {
    expected => 0,
    coords   => [ 11, 11, 10, 30, 0 ],
    comment  => "just 5' of 2,7 inside 3'"
  },
  {
    expected => 0,
    coords   => [ 18, 18, 10, 30, 0 ],
    comment  => "just 3' of 2,7 inside 3'"
  },
  
  # 7,2 inside 5'
  {
    expected => 1,
    coords   => [ 23, 23, 10, 30, 0 ],
    comment  => "just within 7,2 inside 5'"
  },
  {
    expected => 1,
    coords   => [ 24, 26, 10, 30, 0 ],
    comment  => "compeletely within 7,2 inside 5'"
  },
  {
    expected => 1,
    coords   => [ 22, 24, 10, 30, 0 ],
    comment  => "part 5' overlap 7,2 inside 5'"
  },
  {
    expected => 1,
    coords   => [ 26, 29, 10, 30, 0 ],
    comment  => "part 3' overlap 7,2 inside 5'"
  },
  {
    expected => 0,
    coords   => [ 22, 22, 10, 30, 0 ],
    comment  => "just 5' of 7,2 inside 5'"
  },
  {
    expected => 0,
    coords   => [ 29, 29, 10, 30, 0 ],
    comment  => "just 3' of 7,2 inside 5'"
  },
  
  # insertion
  {
    expected => 1,
    coords   => [ 10,  9, 10, 30, 1 ],
    comment  => "insertion intron start"
  },
  {
    expected => 1,
    coords   => [ 12, 11, 10, 30, 1 ],
    comment  => "insertion intron start + 2"
  },
  {
    expected => 1,
    coords   => [  9,  8, 10, 30, 1 ],
    comment  => "insertion intron start 5'"
  },
  {
    expected => 0,
    coords   => [ 11, 10, 10, 30, 1 ],
    comment  => "insertion between intron start, intron start + 2"
  },
  {
    expected => 1,
    coords   => [ 13, 12, 10, 30, 1 ],
    comment  => "insertion intron start + 2 3'"
  },
  {
    expected => 1,
    coords   => [ 31, 30, 10, 30, 1 ],
    comment  => "insertion intron end"
  },
  {
    expected => 1,
    coords   => [ 33, 32, 10, 30, 1 ],
    comment  => "insertion intron end + 2"
  },
  {
    expected => 0,
    coords   => [ 30, 29, 10, 30, 1 ],
    comment  => "insertion intron end 5'"
  },
  {
    expected => 1,
    coords   => [ 32, 31, 10, 30, 1 ],
    comment  => "insertion between intron end, intron end + 2"
  },
  {
    expected => 0,
    coords   => [ 34, 33, 10, 30, 1 ],
    comment  => "insertion intron end + 2 3'"
  },
);

foreach my $c(@coords) {
  # my $p_res   = Bio::EnsEMBL::Variation::Utils::VariationEffect::_intron_overlap_perl(@{$c->{coords}}) || 0;
  my $c_res   = Bio::EnsEMBL::Variation::Utils::VariationEffect::_intron_overlap(@{$c->{coords}}) || 0;
  my $exp     = $c->{expected};
  my $comment = $c->{comment} || " coords: ".join(", ", @{$c->{coords}});
  
  # is($p_res, $exp, "perl intron_overlap $comment");
  is($c_res, $exp, "C    intron_overlap $comment");
}


## TEST DB STUFF

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');
my $vdba = $multi->get_DBAdaptor('variation');
my $cdba = $multi->get_DBAdaptor('core');


my $ta = $cdba->get_TranscriptAdaptor;

my $transcript_tests;

# check a forward strand coding transcript

my $tf = $ta->fetch_by_stable_id('ENST00000360027');

my $t_start = $tf->seq_region_start;
my $t_end   = $tf->seq_region_end;

my $cds_start = $tf->coding_region_start;
my $cds_end   = $tf->coding_region_end;

my $first_intron = $tf->get_all_Introns->[0];

my $intron_start = $first_intron->seq_region_start;
my $intron_end   = $first_intron->seq_region_end;

my $second_exon = $tf->get_all_Exons->[1];

my $exon_start = $second_exon->seq_region_start;
my $exon_end   = $second_exon->seq_region_end;

$transcript_tests->{$tf->stable_id}->{transcript} = $tf;

$transcript_tests->{$tf->stable_id}->{tests} = [
        
    # check the boundaries of the upstream and downstream calls
    
    {
        start   => $t_start - 5001,
        end     => $t_start - 5001,
        effects => [ qw(intergenic_variant) ],
    }, {
        start   => $t_start - 5000,
        end     => $t_start - 5000,
        effects => [ qw(upstream_gene_variant) ],
    }, {
        start   => $t_start - 2001,
        end     => $t_start - 2001,
        effects => [ qw(upstream_gene_variant) ],
    }, {
        start   => $t_start - 2000,
        end     => $t_start - 2000,
        effects => [ qw(upstream_gene_variant) ],
    },{
        start   => $t_start - 1,
        end     => $t_start - 1,
        effects => [ qw(upstream_gene_variant) ],
    }, {
        comment => 'an insertion just before the start is upstream',
        alleles => 'A',
        start   => $t_start,
        end     => $t_start - 1,
        effects => [ qw(upstream_gene_variant) ],
    }, {
        comment => 'an insertion just after the end is downstream',
        alleles => 'A',
        start   => $t_end+1,
        end     => $t_end,
        effects => [ qw(downstream_gene_variant) ],
    }, {
        start   => $t_end + 1,
        end     => $t_end + 1,
        effects => [ qw(downstream_gene_variant) ],
    }, {
        start   => $t_end + 500,
        end     => $t_end + 500,
        effects => [ qw(downstream_gene_variant) ],
    }, {
        start   => $t_end + 501,
        end     => $t_end + 501,
        effects => [ qw(downstream_gene_variant) ],
    }, {   
        start   => $t_end + 5000,
        end     => $t_end + 5000,
        effects => [ qw(downstream_gene_variant) ],
    }, {   
        start   => $t_end + 5001,
        end     => $t_end + 5001,
        effects => [ qw(intergenic_variant) ],
    },

    # check the UTR calls
    
    {
        start   => $t_start,
        end     => $t_start,
        effects => [qw(5_prime_UTR_variant)],
    }, {
        comment => 'an insertion between the first 2 bases is UTR',
        alleles => 'A',
        start   => $t_start + 1,
        end     => $t_start,
        effects => [ qw(5_prime_UTR_variant) ],
    }, {
        start   => $cds_start-1,
        end     => $cds_start-1,
        effects => [qw(5_prime_UTR_variant)],
    }, {
        comment => 'an insertion just before the cds start is UTR',
        alleles => 'T',
        start   => $cds_start, 
        end     => $cds_start-1,
        effects => [qw(5_prime_UTR_variant)],
    }, {
        comment => 'an insertion just after the cds end is UTR',
        alleles => 'A',
        start   => $cds_end+1, 
        end     => $cds_end,
        effects => [qw(3_prime_UTR_variant)],
    }, {
        start   => $cds_end+1,
        end     => $cds_end+1,
        effects => [qw(3_prime_UTR_variant)],
    }, {
        start   => $t_end,
        end     => $t_end,
        effects => [qw(3_prime_UTR_variant)],
    },

    # check the introns & splice sites
    
    {
        start   => $intron_start-4,
        end     => $intron_start-4,
        effects => [qw(missense_variant)],
    }, {
        start   => $intron_start-3,
        end     => $intron_start-3,
        effects => [qw(splice_region_variant synonymous_variant)],
    }, {
        start   => $intron_start,
        end     => $intron_start,
        effects => [qw(splice_donor_variant)],
    }, {
        start   => $intron_start+1,
        end     => $intron_start+1,
        effects => [qw(splice_donor_variant)],
    }, {
        start   => $intron_start+2,
        end     => $intron_start+2,
        effects => [qw(splice_donor_region_variant intron_variant)],
    }, {
        start   => $intron_start+7,
        end     => $intron_start+7,
        effects => [qw(splice_region_variant intron_variant)],
    }, {
        start   => $intron_start+8,
        end     => $intron_start+8,
        effects => [qw(intron_variant)],
    }, {
        comment => 'an insertion between the last exon base and the first intron base is not essential',
        alleles => 'A',
        start   => $intron_start,
        end     => $intron_start-1,
        effects => [qw(splice_region_variant frameshift_variant)],
    }, {
        comment => 'an insertion between the first two bases of an intron is in the donor',
        alleles => 'A',
        start   => $intron_start+1,
        end     => $intron_start,
        effects => [qw(splice_donor_variant)],
    }, {
        comment => 'insertion between bases 2 & 3 of an intron is splice_region',
        alleles => 'A',
        start   => $intron_start+2,
        end     => $intron_start+1,
        effects => [qw(splice_region_variant intron_variant)],
    }, {
        comment => 'insertion between bases 7 & 8 is still splice_region',
        alleles => 'A',
        start   => $intron_start+7,
        end     => $intron_start+6,
        effects => [qw(splice_region_variant intron_variant)],
    }, {
        comment => 'insertion between bases 8 & 9 is just an intron_variant',
        alleles => 'A',
        start   => $intron_start+8,
        end     => $intron_start+7,
        effects => [qw(intron_variant)],
    }, {
        start   => $intron_end - 8,
        end     => $intron_end - 8,
        effects => [qw(intron_variant splice_polypyrimidine_tract_variant)],
    }, {
        start   => $intron_end - 7,
        end     => $intron_end - 7,
        effects => [qw(splice_region_variant intron_variant splice_polypyrimidine_tract_variant)],
    }, {
        start   => $intron_end - 2,
        end     => $intron_end - 2,
        effects => [qw(splice_region_variant intron_variant splice_polypyrimidine_tract_variant)],
    }, {
        start   => $intron_end - 1,
        end     => $intron_end - 1,
        effects => [qw(splice_acceptor_variant)],
    }, {
        start   => $intron_end,
        end     => $intron_end,
        effects => [qw(splice_acceptor_variant)],
    }, {
        start   => $intron_end+1,
        end     => $intron_end+1,
        effects => [qw(splice_region_variant synonymous_variant)],
    }, {
        start   => $intron_end+3,
        end     => $intron_end+3,
        effects => [qw(splice_region_variant missense_variant)],
    }, {
        start   => $intron_end+4,
        end     => $intron_end+4,
        effects => [qw(stop_gained)],
    }, {
        comment => 'an insertion between the last intron base and the first exon base is not essential',
        alleles => 'A',
        start   => $intron_end+1,
        end     => $intron_end,
        effects => [qw(splice_region_variant frameshift_variant)],
    }, {
        comment => 'an insertion between the last two bases of an intron is in the acceptor',
        alleles => 'A',
        start   => $intron_end,
        end     => $intron_end-1,
        effects => [qw(splice_acceptor_variant)],
    }, {
        comment => 'insertion between last bases 2 & 3 of an intron is splice_region',
        alleles => 'T',
        start   => $intron_end-1,
        end     => $intron_end-2,
        effects => [qw(splice_region_variant splice_polypyrimidine_tract_variant intron_variant)],
    }, {
        comment => 'insertion between last bases 7 & 8 is still splice_region',
        alleles => 'A',
        start   => $intron_end-6,
        end     => $intron_end-7,
        effects => [qw(splice_region_variant intron_variant splice_polypyrimidine_tract_variant)],
    }, {
        comment => 'insertion between last bases 8 & 9 is just an intron_variant and splice_polypyrimidine_tract_variant',
        alleles => 'A',
        start   => $intron_end-7,
        end     => $intron_end-8,
        effects => [qw(intron_variant splice_polypyrimidine_tract_variant)],
    }, {
        comment => 'whole exon deletion',
        alleles => '-',
        start   => $exon_start-10,
        end     => $exon_end+10,
        effects => [qw(intron_variant splice_donor_5th_base_variant splice_acceptor_variant splice_donor_variant coding_sequence_variant )],
    }, {
        comment => 'long sequence var where middle is identical but overlaps splice site',
        alleles => 'ATGTACTGCCTATGTGTGCTGTGAGTATGATACGGTGGACT',
        start   => $intron_start - 20,
        end     => $intron_start + 20,
        effects => [qw(intron_variant coding_sequence_variant)],
    },

    # check the CDS 

    {
        alleles => 'G',
        start   => $cds_start,
        end     => $cds_start,
        effects => [qw(start_lost)],
    }, {
        alleles => 'G',
        start   => $cds_start+1,
        end     => $cds_start+1,
        effects => [qw(start_lost)],
    }, {
        alleles => 'C',
        start   => $cds_start+2,
        end     => $cds_start+2,
        effects => [qw(start_lost)],
    },  {
        alleles => 'C',
        start   => $cds_start+3,
        end     => $cds_start+3,
        effects => [qw(missense_variant)],
 },  {
        alleles => 'CGGTGT',
        start   => $cds_start+3,
        end     => $cds_start+5,
        effects => [qw(protein_altering_variant)],
    }, {
        alleles => 'CCC',
        start   => $cds_start+3,
        end     => $cds_start+2,
        effects => [qw(inframe_insertion)],
    }, {
        alleles => 'CCC',
        start   => $cds_start+2,
        end     => $cds_start+1,
        effects => [qw(start_lost)],
    }, {
        alleles => 'AGG',
        start   => $cds_start+2,
        end     => $cds_start+1,
        effects => [qw(start_lost)],
    }, {
        alleles => 'CAT',
        start   => $cds_start+2,
        end     => $cds_start+1,
        effects => [qw(start_lost start_retained_variant)],
    }, {
        alleles => 'GCA',
        start   => $cds_start+2,
        end     => $cds_start+1,
        effects => [qw(inframe_insertion start_retained_variant)],
    }, {
        alleles => '-',
        start   => $cds_start+3,
        end     => $cds_start+5,
        effects => [qw(inframe_deletion)],
        pep_alleles => 'D/-',
    }, {
        alleles => '-',
        start   => $cds_start+4,
        end     => $cds_start+6,
        effects => [qw(inframe_deletion)],
        pep_alleles => 'DA/A',
    }, {
        alleles => 'GAT',
        start   => $cds_start+3,
        end     => $cds_start+5,
        effects => [qw(synonymous_variant)],
        pep_alleles => 'D/D',
    }, {
        alleles => 'GATACA',
        start   => $cds_start+3,
        end     => $cds_start+8,
        effects => [qw(missense_variant)],
        pep_alleles => 'DA/DT',
    }, {
        alleles => 'CAT',
        start   => $cds_start-1,
        end     => $cds_start+2,
        effects => [qw(5_prime_UTR_variant start_retained_variant)],
        pep_alleles => '',
    }, {
        alleles => 'G',
        start   => $cds_start+4,
        end     => $cds_start+3,
        effects => [qw(frameshift_variant)],
    }, {
        alleles => 'GT',
        start   => $cds_start+4,
        end     => $cds_start+3,
        effects => [qw(frameshift_variant)],
    }, {
        alleles => 'GTAG',
        start   => $cds_start+4,
        end     => $cds_start+3,
        effects => [qw(frameshift_variant)],
    }, {
        alleles => '-',
        start   => $cds_start+3,
        end     => $cds_start+3,
        effects => [qw(frameshift_variant)],
    }, {
        alleles => '-',
        start   => $cds_start+3,
        end     => $cds_start+4,
        effects => [qw(frameshift_variant)],
    }, {
        alleles => '-',
        start   => $cds_start+3,
        end     => $cds_start+6,
        effects => [qw(frameshift_variant)],
    }, {
        alleles => 'G',
        start   => $cds_end-2,
        end     => $cds_end-2,
        effects => [qw(stop_lost)],
    }, {
        alleles => 'A',
        start   => $cds_end-1,
        end     => $cds_end-1,
        effects => [qw(stop_retained_variant)],
    }, {
        alleles => 'C',
        start   => $cds_end,
        end     => $cds_end,
        effects => [qw(stop_lost)],
    }, {
        alleles => 'AAG',
        start   => $cds_end-1,
        end     => $cds_end-2,
        effects => [qw(inframe_insertion stop_retained_variant )],
    }, {
        alleles => '-',
        start   => $cds_end-2,
        end     => $cds_end,
        effects => [qw(stop_lost inframe_deletion)],
    }, {
        alleles => 'TAA',
        start   => $cds_end-2,
        end     => $cds_end,
        effects => [qw(stop_retained_variant)],
    }, {
        alleles => 'GGG',
        start   => $cds_end-2,
        end     => $cds_end,
        effects => [qw(stop_lost)],
    }, {
        comment => 'insertion in stop codon that shifts stop downstream and has additional bases beyond stop',
        alleles => 'CTGAGG',
        start   => $cds_end,
        end     => $cds_end-1,
        effects => [qw(inframe_insertion)],
    },
    
    # =============================================================================
    # TEST CASES FOR GITHUB ISSUE ensembl-vep#1710
    # =============================================================================
    # Bug: VEP 112+ incorrectly predicts stop_retained_variant instead of stop_gained
    # for insertions that introduce a stop codon when the reference has no stop.
    #
    # Root cause: In ref_eq_alt_sequence(), condition 1 was:
    #   ($ref_pep eq substr($alt_pep, 0, 1) && $alt_pep =~ /\*/)
    # This returned stop_retained when first AA matched and alt had stop, WITHOUT
    # checking if ref also had a stop codon.
    #
    # Fix: Added check that ref must also have stop codon for stop_retained.
    #
    # The following tests cover:
    # 1. Issue #1710 main case: insertion with embedded stop -> should be stop_gained
    # 2. Edge cases: various combinations of ref/alt with/without stop codons
    # 3. Regression tests: existing stop_retained cases should still work
    # =============================================================================
    
    # ---------------------------------------------------------------------------
    # ISSUE #1710 REGRESSION TESTS: Insertions with embedded stop codons
    # These are the PRIMARY bug cases that Issue #1710 reported
    # ---------------------------------------------------------------------------
    {
        # Issue #1710 Case 1: Inframe insertion that introduces a NEW stop codon
        # Example from issue: 3:56591278-56591278 T>TGGGGTAAGCA (L -> LG*AX)
        # ref_pep = single amino acid (no stop), alt_pep = amino acids with embedded stop
        # This should be stop_gained, NOT stop_retained
        # The ref has NO stop codon, so there is nothing to "retain"
        #
        # COORDINATES EXPLANATION:
        # - Insert BETWEEN $cds_end-3 (last base of penultimate codon) and $cds_end-2 (first base of stop)
        # - This places the insertion at a codon boundary so TAAGGG becomes codons TAA + GGG
        # - Original: ...XXX|TAA (where XXX is penultimate codon, | is insertion point, TAA is stop)
        # - After: ...XXX|TAAGGG|TAA -> codons: XXX TAA GGG TAA
        # - The inserted TAA becomes a new in-frame stop codon = stop_gained
        comment => 'GitHub Issue #1710: inframe insertion with embedded stop should be stop_gained not stop_retained',
        alleles => 'TAAGGG',  # Inserting TAA (stop) + GGG (Gly) - TAA first so it's in-frame as stop
        start   => $cds_end-2,  # First base of original stop codon
        end     => $cds_end-3,  # Last base of penultimate codon (insertion between them)
        effects => [qw(stop_gained inframe_insertion)],
    }, {
        # Issue #1710 Case 2: Multi-codon insertion with stop in middle
        # Insert at codon boundary so the embedded stop is in-frame
        # ref has no stop at this position, alt gains a stop -> stop_gained
        #
        # COORDINATES: Insert between $cds_end-6 and $cds_end-5 (codon boundary)
        # This is a 9bp (3 codon) insertion, so still inframe
        comment => 'GitHub Issue #1710: multi-codon insertion with stop in middle should be stop_gained',
        alleles => 'GGGTAGGGG',  # GGG (Gly) + TAG (stop) + GGG (Gly - but after stop, won't translate)
        start   => $cds_end-5,  # First base of second-to-last coding codon
        end     => $cds_end-6,  # Last base of third-to-last coding codon (insertion between them)
        effects => [qw(stop_gained inframe_insertion)],
    },
    
    # ---------------------------------------------------------------------------
    # EDGE CASE: Single nucleotide changes at stop codon positions
    # These should continue to work correctly after the fix
    # ---------------------------------------------------------------------------
    {
        # Substitution changing last base of stop codon
        # The stop codon is TGA. Changing 'A' to 'G' creates TGG (Trp) = stop_lost
        # NOTE: This test was previously incorrectly expecting stop_retained,
        # assuming the stop was TAA and this was a synonymous change.
        # The actual stop codon in the test transcript is TGA.
        comment => 'Edge case: TGA last base change to G creates TGG = stop_lost',
        alleles => 'G',
        start   => $cds_end,  # Last base of stop codon
        end     => $cds_end,
        effects => [qw(stop_lost)],
    },
    
    # ---------------------------------------------------------------------------
    # EDGE CASE: Deletions affecting stop codons
    # ---------------------------------------------------------------------------
    {
        # Deletion that removes stop codon entirely
        # ref has stop, alt has no stop -> stop_lost
        comment => 'Edge case: deletion removing entire stop codon should be stop_lost',
        alleles => '-',
        start   => $cds_end-2,
        end     => $cds_end,
        effects => [qw(stop_lost inframe_deletion)],
    },
    
    # ---------------------------------------------------------------------------
    # EDGE CASE: Complex insertions near stop codon
    # ---------------------------------------------------------------------------
    {
        # Insertion right before stop that does NOT introduce a new stop
        # Insert at codon boundary (between penultimate and stop codon)
        # 3 bases = inframe insertion, GGG = Gly (not a stop)
        # ref has stop, insertion pushes stop 3bp downstream -> still stop_retained? No, this is inframe_insertion only
        comment => 'Edge case: insertion before stop without new stop should not be stop_retained',
        alleles => 'GGG',  # Three bases (Gly) without stop
        start   => $cds_end-2,  # First base of stop codon
        end     => $cds_end-3,  # Last base of penultimate codon (insertion between them)
        effects => [qw(inframe_insertion)],
    }, {
        # Insertion of TAA within the stop codon
        # This inserts TAA between the 2nd and 3rd base of the stop codon TGA
        # Original: TGA (stop) -> After insertion: TGTAAA (inserting TAA between G and A)
        # But the CI shows: tga/tgTAAa which translates to */CK
        # The inserted TAA disrupts the original stop codon frame
        # Result: ref has stop (*), alt has CK (no stop) = stop_lost
        # VEP behavior: Only stop_lost is returned (not also inframe_insertion).
        # When a stop codon is lost, VEP doesn't additionally report inframe_insertion.
        comment => 'Edge case: TAA insertion within stop codon disrupts stop = stop_lost',
        alleles => 'TAA',  # Inserting 3 bases within stop codon
        start   => $cds_end,
        end     => $cds_end-1,
        effects => [qw(stop_lost)],
    },
    
    # ---------------------------------------------------------------------------
    # EDGE CASE: Substitutions involving stop codons
    # These test the boundary between stop_gained, stop_lost, and stop_retained
    # ---------------------------------------------------------------------------
    {
        # Substitution changing stop codon to another stop codon (TAA -> TAG)
        # Both ref and alt have stop -> stop_retained
        comment => 'Edge case: stop codon to different stop codon should be stop_retained',
        alleles => 'TAG',
        start   => $cds_end-2,
        end     => $cds_end,
        effects => [qw(stop_retained_variant)],
    }, {
        # Substitution changing stop codon to non-stop (TAA -> GAA)
        # ref has stop, alt has no stop -> stop_lost
        comment => 'Edge case: stop codon to non-stop should be stop_lost',
        alleles => 'GAA',
        start   => $cds_end-2,
        end     => $cds_end,
        effects => [qw(stop_lost)],
    }, {
        # Substitution in coding region that creates a new stop (creates premature stop)
        # ref has no stop at this position, alt has stop -> stop_gained
        comment => 'Edge case: non-stop to stop codon should be stop_gained',
        alleles => 'TAA',
        start   => $cds_end-5,  # Before the actual stop codon
        end     => $cds_end-3,
        effects => [qw(stop_gained)],
    },
    
    # ---------------------------------------------------------------------------
    # BIOLOGICALLY UNUSUAL CASES (may not occur naturally but test logic)
    # ---------------------------------------------------------------------------
    {
        # Very long insertion with stop near the beginning
        # Tests that we correctly identify stop_gained even with trailing sequence
        # 9 bases = inframe (3 codons), insert at codon boundary
        # Codons: TAA (stop) + GGG (Gly) + AAA (Lys) - stop is first, so stop_gained
        comment => 'Unusual: long insertion with early stop should be stop_gained',
        alleles => 'TAAGGGAAA',  # TAA (stop) + GGG + AAA
        start   => $cds_end-5,  # First base of second-to-last coding codon
        end     => $cds_end-6,  # Last base of third-to-last coding codon (insertion between)
        effects => [qw(stop_gained inframe_insertion)],
    }, {
        # Insertion of just a stop codon (TAA) in coding region
        # 3 bases = inframe, insert at codon boundary
        # ref has no stop here, alt has stop -> stop_gained
        comment => 'Unusual: insertion of bare stop codon should be stop_gained',
        alleles => 'TAA',
        start   => $cds_end-5,  # First base of second-to-last coding codon
        end     => $cds_end-6,  # Last base of third-to-last coding codon (insertion between)
        effects => [qw(stop_gained inframe_insertion)],
    },
    
    # ---------------------------------------------------------------------------
    # NEGATIVE TESTS: Ensure we DON'T incorrectly call stop_retained
    # These are the specific patterns that triggered Issue #1710
    # ---------------------------------------------------------------------------
    {
        # Pattern: L -> LG*AX (ref single AA, alt has embedded stop)
        # This was the exact Issue #1710 bug pattern
        # First char matches (L=L), alt has stop, but ref has NO stop
        # Must NOT be stop_retained, MUST be stop_gained
        # 6 bases = inframe, TGAGGG = TGA (stop) + GGG (Gly)
        # Insert at codon boundary so TGA becomes a proper stop codon
        comment => 'NEGATIVE TEST #1710: single AA to multi-AA-with-stop must NOT be stop_retained',
        alleles => 'TGAGGG',  # TGA (stop) + GGG (Gly) - stop first for in-frame stop_gained
        start   => $cds_end-2,  # First base of stop codon
        end     => $cds_end-3,  # Last base of penultimate codon (insertion between)
        effects => [qw(stop_gained inframe_insertion)],
    },
    
    # =============================================================================
    # END OF ISSUE #1710 TEST CASES
    # =============================================================================

    # =============================================================================
    # ENSVAR-6654 NOTES (from PR #1184)
    # =============================================================================
    # Additional edge cases from PR #1184 (ENSVAR-6654) are incorporated into the
    # code but are difficult to trigger through artificial test coordinates:
    #   - inframe_insertion guard: return 0 when both ref_pep and alt_pep are '*'
    #     (e.g., codons TAA/TAAG where alt is 1bp larger but still just a stop)
    #   - stop_lost/stop_retained: fall through to sequence-level analysis when
    #     alt_pep contains 'X' (unknown/incomplete amino acid translation)
    #   - _overlaps_stop_codon: insertion coordinate fix (extends cdna_end by
    #     insertion length to detect overlap correctly)
    #   - ref_eq_alt_sequence: condition improvements (remove redundant condition,
    #     check trailing stop semantically, simplify index comparison)
    # These changes were validated by the Ensembl team against the full GRCh38
    # variant database with no difference in variant consequences.
    # =============================================================================

    # =============================================================================
    # TEST CASES FOR GITHUB ISSUE ensembl-vep#1710 - BUG 2 (FORWARD STRAND)
    # =============================================================================
    # Bug 2: Frameshift insertions/deletions at the stop codon incorrectly return
    # stop_retained_variant instead of frameshift_variant + stop_lost.
    #
    # Root cause: _ins_del_stop_altered() didn't detect frameshifts - it only
    # checked if the codon at the original stop position was still a stop,
    # which is semantically wrong for frameshifts where the reading frame shifts.
    #
    # Fix: Added frameshift detection (non-3n length change) in _ins_del_stop_altered()
    # to return TRUE (stop IS altered) for all frameshift variants at stop codon.
    #
    # Stop codon positions (forward strand): $cds_end-2, $cds_end-1, $cds_end
    # Insertion notation: start => X, end => X-1 means insert BETWEEN X-1 and X
    # =============================================================================

    # ---------------------------------------------------------------------------
    # POSITION-BASED INSERTION TESTS
    # ---------------------------------------------------------------------------
    {
        # 1bp insertion BEFORE stop codon (in last coding codon, not stop itself)
        # This causes frameshift but doesn't overlap stop, so no stop_lost
        comment => 'Bug2: 1bp insertion BEFORE stop codon (frameshift in coding)',
        alleles => 'A',
        start   => $cds_end-2,
        end     => $cds_end-3,
        effects => [qw(frameshift_variant)],
    }, {
        # 1bp insertion BETWEEN 1st and 2nd base of stop (overlaps stop)
        comment => 'Bug2: 1bp insertion between 1st-2nd base of stop = frameshift+stop_lost',
        alleles => 'A',
        start   => $cds_end-1,
        end     => $cds_end-2,
        effects => [qw(frameshift_variant stop_lost)],
    }, {
        # 1bp insertion BETWEEN 2nd and 3rd base of stop (overlaps stop)
        comment => 'Bug2: 1bp insertion between 2nd-3rd base of stop = frameshift+stop_lost',
        alleles => 'A',
        start   => $cds_end,
        end     => $cds_end-1,
        effects => [qw(frameshift_variant stop_lost)],
    }, {
        # 1bp insertion AFTER stop codon (at CDS/UTR boundary)
        # This insertion is in the 3' UTR, right after the stop codon ends.
        # Since it's in the UTR (not the CDS), it doesn't cause a frameshift
        # or affect the stop codon - the stop has already terminated translation.
        # NOTE: Original expectation of frameshift+stop_lost was incorrect.
        # The variant is purely in the UTR.
        comment => 'Bug2: 1bp insertion at CDS/UTR boundary = UTR variant only',
        alleles => 'A',
        start   => $cds_end+1,
        end     => $cds_end,
        effects => [qw(3_prime_UTR_variant)],
    }, {
        # 1bp insertion purely in UTR (no stop codon overlap)
        comment => 'Bug2: 1bp insertion in 3-prime UTR only',
        alleles => 'A',
        start   => $cds_end+2,
        end     => $cds_end+1,
        effects => [qw(3_prime_UTR_variant)],
    },

    # ---------------------------------------------------------------------------
    # SIZE-BASED INSERTION TESTS (at 2nd-3rd base position)
    # ---------------------------------------------------------------------------
    {
        # 2bp insertion = frameshift (2 % 3 != 0)
        comment => 'Bug2: 2bp insertion at stop = frameshift+stop_lost',
        alleles => 'AT',
        start   => $cds_end,
        end     => $cds_end-1,
        effects => [qw(frameshift_variant stop_lost)],
    }, {
        # 4bp insertion = frameshift (4 % 3 != 0)
        comment => 'Bug2: 4bp insertion at stop = frameshift+stop_lost',
        alleles => 'ATCG',
        start   => $cds_end,
        end     => $cds_end-1,
        effects => [qw(frameshift_variant stop_lost)],
    }, {
        # 5bp insertion = frameshift (5 % 3 != 0)
        comment => 'Bug2: 5bp insertion at stop = frameshift+stop_lost',
        alleles => 'ATCGA',
        start   => $cds_end,
        end     => $cds_end-1,
        effects => [qw(frameshift_variant stop_lost)],
    }, {
        # 7bp insertion = frameshift (7 % 3 != 0)
        comment => 'Bug2: 7bp insertion at stop = frameshift+stop_lost',
        alleles => 'ATCGATC',
        start   => $cds_end,
        end     => $cds_end-1,
        effects => [qw(frameshift_variant stop_lost)],
    },

    # ---------------------------------------------------------------------------
    # DELETION TESTS (frameshift deletions at stop codon)
    # ---------------------------------------------------------------------------
    {
        # 1bp deletion of 1st base of stop
        comment => 'Bug2: 1bp deletion of 1st base of stop = frameshift+stop_lost',
        alleles => '-',
        start   => $cds_end-2,
        end     => $cds_end-2,
        effects => [qw(frameshift_variant stop_lost)],
    }, {
        # 1bp deletion of 2nd base of stop
        comment => 'Bug2: 1bp deletion of 2nd base of stop = frameshift+stop_lost',
        alleles => '-',
        start   => $cds_end-1,
        end     => $cds_end-1,
        effects => [qw(frameshift_variant stop_lost)],
    }, {
        # 1bp deletion of 3rd base of stop
        comment => 'Bug2: 1bp deletion of 3rd base of stop = frameshift+stop_lost',
        alleles => '-',
        start   => $cds_end,
        end     => $cds_end,
        effects => [qw(frameshift_variant stop_lost)],
    }, {
        # 2bp deletion from stop (1st+2nd base)
        comment => 'Bug2: 2bp deletion from stop = frameshift+stop_lost',
        alleles => '-',
        start   => $cds_end-2,
        end     => $cds_end-1,
        effects => [qw(frameshift_variant stop_lost)],
    }, {
        # 4bp deletion spanning stop into UTR
        # Deletes 2 bases of stop codon + 2 bases of UTR
        # VEP behavior: stop_lost + 3_prime_UTR_variant (no frameshift reported
        # for deletions that span CDS/UTR boundary into UTR)
        comment => 'Bug2: 4bp deletion spanning stop into UTR = stop_lost+UTR',
        alleles => '-',
        start   => $cds_end-1,
        end     => $cds_end+2,
        effects => [qw(3_prime_UTR_variant stop_lost)],
    }, {
        # 4bp deletion from coding into stop
        comment => 'Bug2: 4bp deletion from coding into stop = frameshift+stop_lost',
        alleles => '-',
        start   => $cds_end-4,
        end     => $cds_end-1,
        effects => [qw(frameshift_variant stop_lost)],
    },

    # ---------------------------------------------------------------------------
    # NEGATIVE TESTS: In-frame variants (should NOT trigger Bug2 fix)
    # These verify we didn't break existing in-frame handling
    # ---------------------------------------------------------------------------
    {
        # 3bp insertion = in-frame, uses existing stop codon check logic
        comment => 'Bug2 NEGATIVE: 3bp insertion at stop = in-frame (existing logic)',
        alleles => 'ATG',
        start   => $cds_end,
        end     => $cds_end-1,
        effects => [qw(inframe_insertion stop_retained_variant)],
    }, {
        # 6bp insertion = in-frame (two codons)
        comment => 'Bug2 NEGATIVE: 6bp insertion at stop = in-frame (existing logic)',
        alleles => 'ATGATG',
        start   => $cds_end,
        end     => $cds_end-1,
        effects => [qw(inframe_insertion stop_retained_variant)],
    }, {
        # 9bp insertion = in-frame (three codons)
        comment => 'Bug2 NEGATIVE: 9bp insertion at stop = in-frame (existing logic)',
        alleles => 'ATGATGATG',
        start   => $cds_end,
        end     => $cds_end-1,
        effects => [qw(inframe_insertion stop_retained_variant)],
    },

    # ---------------------------------------------------------------------------
    # COMPLEX EDGE CASE
    # ---------------------------------------------------------------------------
    {
        # Deletion spanning from coding through stop into UTR (5bp = frameshift)
        # VEP behavior: stop_lost + 3_prime_UTR_variant (frameshift not reported
        # for deletions spanning CDS/UTR boundary)
        comment => 'Bug2 Complex: 5bp deletion from coding through stop into UTR',
        alleles => '-',
        start   => $cds_end-3,
        end     => $cds_end+1,
        effects => [qw(3_prime_UTR_variant stop_lost)],
    },

    # =============================================================================
    # END OF ISSUE #1710 BUG 2 FORWARD STRAND TEST CASES
    # =============================================================================

    {
        comment => 'a wierd allele string',
        alleles => 'HGMD_MUTATION',
        start   => $cds_end-10,
        end     => $cds_end-11,
        effects => [qw(coding_sequence_variant)],
    }, {
        comment => 'an ambiguous allele string',
        alleles => 'W',
        start   => $cds_end-10,
        end     => $cds_end-10,
        effects => [qw(coding_sequence_variant)],
    }, {
        comment => 'an ambiguous insertion',
        alleles => 'W',
        start   => $cds_end-9,
        end     => $cds_end-10,
        effects => [qw(frameshift_variant)],
    }, {
        comment => 'a specified length insertion',
        alleles => '2 BP INSERTION',
        start   => $cds_end-9,
        end     => $cds_end-10,
        effects => [qw(frameshift_variant)],
    }, {
        comment => 'an inframe specified length insertion',
        alleles => '3 BP INSERTION',
        start   => $cds_end-9,
        end     => $cds_end-10,
        effects => [qw(coding_sequence_variant)],
    }, {
        comment => 'delete the last codon of an exon - shifting into splice donor region',
        alleles => '-',
        start   => $intron_start-3,
        end     => $intron_start-1,
        effects => [qw(coding_sequence_variant splice_donor_variant)],
        no_shift => 0,
    }, 
    

    # check the complex calls
   
    {
        alleles => '-',
        start   => $intron_start-3,
        end     => $intron_start+2,
        effects => [qw( splice_donor_variant splice_donor_region_variant coding_sequence_variant intron_variant)],
    }, {
        alleles => '-',
        start   => $intron_end-2,
        end     => $intron_end+3,
        effects => [qw( splice_acceptor_variant coding_sequence_variant intron_variant)],
    }, {
        comment => 'deletion overlapping UTR and start site, start lost 1',
        alleles => '-',
        start   => $cds_start-3,
        end     => $cds_start+2,
        effects => [qw( 5_prime_UTR_variant start_lost)],
    },  {
        comment => 'deletion overlapping UTR and start site, start lost 2',
        alleles => '-',
        start   => $cds_start-1, 
        end     => $cds_start+1,
        effects => [qw(5_prime_UTR_variant start_lost)],
    }, {
        comment => 'unbalanced sub overlapping UTR and start site, start retained',
        alleles => '-',
        start   => $cds_start-4,
        end     => $cds_start,
        effects => [qw(5_prime_UTR_variant start_retained_variant start_lost)],
    }, {
        comment => 'deletion overlapping STOP and 3\' UTR, stop retained - shifted into solely 3\' UTR',
        alleles => '-',
        start   => $cds_end-1,
        end     => $cds_end+1,
        no_shift => 0,
        effects => [qw( 3_prime_UTR_variant)],
    }, {
        # Deletion: 4 bases from $cds_end-1 to $cds_end+2 (2 of stop + 2 of UTR)
        # Original expectation was stop_retained, but VEP returns stop_lost.
        # This makes sense: deleting 2 bases of the 3-base stop codon changes the reading
        # frame and typically results in a different codon (not a stop).
        comment => 'deletion overlapping STOP and 3\' UTR, stop lost (not retained)',
        alleles => '-',
        start   => $cds_end-1,
        end     => $cds_end+2,
        effects => [qw( 3_prime_UTR_variant stop_lost)],
    }, {
        comment => 'deletion overlapping STOP and 3\' UTR, stop lost',
        alleles => 'C',
        start   => $cds_end-1,
        end     => $cds_end+2,
        effects => [qw( 3_prime_UTR_variant stop_lost)],
    }, {
        comment => 'deletion overlapping STOP and beyond 3\' UTR, stop lost',
        alleles => '-',
        start   => $cds_end-1,
        end     => $cds_end+97,
        effects => [qw( 3_prime_UTR_variant stop_lost)],
    },  

];

####################################################################################

# now do the same for a reverse strand transcript

my $tr = $ta->fetch_by_stable_id('ENST00000368312');

$transcript_tests->{$tr->stable_id}->{transcript} = $tr;

$t_start = $tr->seq_region_start;
$t_end   = $tr->seq_region_end;

$cds_start = $tr->coding_region_start;
$cds_end   = $tr->coding_region_end;

$first_intron = $tr->get_all_Introns->[0];

$intron_start = $first_intron->seq_region_start;
$intron_end   = $first_intron->seq_region_end;

$transcript_tests->{$tr->stable_id}->{tests} = [
        
    # check the boundaries of the upstream and downstream calls
    
    {
        start   => $t_end + 5001,
        end     => $t_end + 5001,
        effects => [ qw(intergenic_variant) ],
    }, {
        start   => $t_end + 5000,
        end     => $t_end + 5000,
        effects => [ qw(upstream_gene_variant) ],
    }, {
        start   => $t_end + 2001,
        end     => $t_end + 2001,
        effects => [ qw(upstream_gene_variant) ],
    }, {
        start   => $t_end + 2000,
        end     => $t_end + 2000,
        effects => [ qw(upstream_gene_variant) ],
    },{
        start   => $t_end + 1,
        end     => $t_end + 1,
        effects => [ qw(upstream_gene_variant) ],
    }, {
        comment => 'an insertion just before the start is upstream',
        alleles => 'G',
        start   => $t_end + 1,
        end     => $t_end,
        effects => [ qw(upstream_gene_variant) ],
        strand  => -1,
    }, {
        comment => 'an insertion just after the end is downstream',
        alleles => 'A',
        start   => $t_start,
        end     => $t_start - 1,
        effects => [ qw(downstream_gene_variant) ],
    }, {
        start   => $t_start - 1,
        end     => $t_start - 1,
        effects => [ qw(downstream_gene_variant) ],
    }, {
        start   => $t_start - 500,
        end     => $t_start - 500,
        effects => [ qw(downstream_gene_variant) ],
    }, {
        start   => $t_start - 501,
        end     => $t_start - 501,
        effects => [ qw(downstream_gene_variant) ],
    }, {   
        start   => $t_start - 5000,
        end     => $t_start - 5000,
        effects => [ qw(downstream_gene_variant) ],
    }, {   
        start   => $t_start - 5001,
        end     => $t_start - 5001,
        effects => [ qw(intergenic_variant) ],
    },

    # check the UTR calls
    
    {
        start   => $t_end,
        end     => $t_end,
        effects => [qw(5_prime_UTR_variant)],
    }, {
        comment => 'an insertion between the first 2 bases is UTR',
        alleles => 'A',
        start   => $t_end,
        end     => $t_end - 1,
        effects => [ qw(5_prime_UTR_variant) ], 
    }, {
        start   => $cds_end + 1,
        end     => $cds_end + 1,
        effects => [qw(5_prime_UTR_variant)],
    }, {
        comment => 'an insertion just before the cds start is UTR',
        alleles => 'A',
        start   => $cds_end + 1, 
        end     => $cds_end,
        effects => [qw(5_prime_UTR_variant)],
    }, {
        comment => 'an insertion just after the cds end is UTR',
        alleles => 'A',
        start   => $cds_start, 
        end     => $cds_start - 1,
        effects => [qw(3_prime_UTR_variant)],
    }, {
        start   => $cds_start - 1,
        end     => $cds_start - 1,
        effects => [qw(3_prime_UTR_variant)],
    }, {
        start   => $t_start,
        end     => $t_start,
        effects => [qw(3_prime_UTR_variant)],
    },

    # check the introns & splice sites
    
    {
        start   => $intron_end + 4,
        end     => $intron_end + 4,
        effects => [qw(synonymous_variant)],
    }, {
        start   => $intron_end + 3,
        end     => $intron_end + 3,
        effects => [qw(splice_region_variant missense_variant)],
    }, {
        start   => $intron_end,
        end     => $intron_end,
        effects => [qw(splice_donor_variant)],
    }, {
        start   => $intron_end - 1,
        end     => $intron_end - 1,
        effects => [qw(splice_donor_variant)],
    }, {
        start   => $intron_end - 2,
        end     => $intron_end - 2,
        effects => [qw(splice_donor_region_variant intron_variant)],
    }, {
        start   => $intron_end - 7,
        end     => $intron_end - 7,
        effects => [qw(splice_region_variant intron_variant)],
    }, {
        start   => $intron_end - 8,
        end     => $intron_end - 8,
        effects => [qw(intron_variant)],
    }, {
        comment => 'an insertion between the last exon base and the first intron base is not essential',
        alleles => 'A',
        start   => $intron_end + 1,
        end     => $intron_end,
        effects => [qw(splice_region_variant frameshift_variant)],
    }, {
        comment => 'an insertion between the first two bases of an intron is in the donor',
        alleles => 'T',
        start   => $intron_end,
        end     => $intron_end - 1,
        effects => [qw(splice_donor_variant)],
    }, {
        comment => 'insertion between bases 2 & 3 of an intron is splice_region',
        alleles => 'T',
        start   => $intron_end - 1,
        end     => $intron_end - 2,
        effects => [qw(splice_region_variant intron_variant)],
    }, {
        comment => 'insertion between bases 7 & 8 is still splice_region',
        alleles => 'A',
        start   => $intron_end - 6,
        end     => $intron_end - 7,
        effects => [qw(splice_region_variant intron_variant)],
    }, {
        comment => 'insertion between bases 8 & 9 is just an intron_variant',
        alleles => 'A',
        start   => $intron_end - 7,
        end     => $intron_end - 8,
        effects => [qw(intron_variant)],
    }, {
        start   => $intron_start + 8,
        end     => $intron_start + 8,
        effects => [qw(intron_variant splice_polypyrimidine_tract_variant)],
    }, {
        start   => $intron_start + 7,
        end     => $intron_start + 7,
        effects => [qw(splice_region_variant intron_variant splice_polypyrimidine_tract_variant)],
    }, {
        start   => $intron_start + 2,
        end     => $intron_start + 2,
        effects => [qw(splice_region_variant intron_variant splice_polypyrimidine_tract_variant)],
    }, {
        start   => $intron_start + 1,
        end     => $intron_start + 1,
        effects => [qw(splice_acceptor_variant)],
    }, {
        start   => $intron_start,
        end     => $intron_start,
        effects => [qw(splice_acceptor_variant)],
    }, {
        start   => $intron_start - 1,
        end     => $intron_start - 1,
        effects => [qw(splice_region_variant missense_variant)],
    }, {
        start   => $intron_start - 3,
        end     => $intron_start - 3,
        effects => [qw(splice_region_variant missense_variant)],
    }, {
        start   => $intron_start - 4,
        end     => $intron_start - 4,
        effects => [qw(missense_variant)],
    }, {
        comment => 'an insertion between the last intron base and the first exon base is not essential',
        alleles => 'A',
        start   => $intron_start,
        end     => $intron_start - 1,
        effects => [qw(splice_region_variant frameshift_variant)],
    }, {
        comment => 'an insertion between the last two bases of an intron is in the acceptor',
        alleles => 'A',
        start   => $intron_start + 1,
        end     => $intron_start,
        effects => [qw(splice_acceptor_variant)],
    }, {
        comment => 'insertion between last bases 2 & 3 of an intron is splice_region',
        alleles => 'A',
        start   => $intron_start + 2,
        end     => $intron_start + 1,
        effects => [qw(splice_region_variant splice_polypyrimidine_tract_variant intron_variant)],
    }, {
        comment => 'insertion between last bases 7 & 8 is still splice_region',
        alleles => 'A',
        start   => $intron_start + 7,
        end     => $intron_start + 6,
        effects => [qw(splice_region_variant intron_variant splice_polypyrimidine_tract_variant)],
    }, {
        comment => 'insertion between last bases 8 & 9 is an intron_variant and splice_polypyrimidine_tract_variant',
        alleles => 'A',
        start   => $intron_start + 8,
        end     => $intron_start + 7,
        effects => [qw(intron_variant splice_polypyrimidine_tract_variant)],
    }, 

    # check the CDS 

    {
        alleles => 'C',
        strand  => -1,
        start   => $cds_end,
        end     => $cds_end,
        effects => [qw(missense_variant)],
    }, {
        alleles => 'G',
        strand  => -1,
        start   => $cds_end - 1,
        end     => $cds_end - 1,
        effects => [qw(missense_variant)],
    }, {
        alleles => 'C',
        strand  => -1,
        start   => $cds_end - 2,
        end     => $cds_end - 2,
        effects => [qw(missense_variant)],
    },  {
        alleles => 'G',
        strand  => -1,
        start   => $cds_end - 3,
        end     => $cds_end - 3,
        effects => [qw(missense_variant)],
    }, {
        alleles => 'GGG',
        strand  => -1,
        start   => $cds_end - 2,
        end     => $cds_end - 3,
        effects => [qw(inframe_insertion)],
    }, {
        alleles => 'GGG',
        strand  => -1,
        start   => $cds_end - 1,
        end     => $cds_end - 2,
        effects => [qw(inframe_insertion)],
    }, {
        alleles => 'AGG',
        strand  => -1,
        start   => $cds_end - 1,
        end     => $cds_end - 2,
        no_shift => 0,
        effects => [qw(protein_altering_variant)],
    }, {
        alleles => '-',
        strand  => -1,
        start   => $cds_end - 5,
        end     => $cds_end - 3,
        effects => [qw(inframe_deletion)],
        pep_alleles => 'L/-',
    }, {
        alleles => 'CTT',
        strand  => -1,
        start   => $cds_end - 5,
        end     => $cds_end - 3,
        effects => [qw(synonymous_variant)],
        pep_alleles => 'L/L',
    }, {
        alleles => 'GATACA',
        strand  => -1,
        start   => $cds_end - 8,
        end     => $cds_end - 3,
        effects => [qw(missense_variant)],
        pep_alleles => 'LT/DT',
    }, {
        alleles => 'G',
        strand  => -1,
        start   => $cds_end - 3,
        end     => $cds_end - 4,
        effects => [qw(frameshift_variant)],
    }, {
        alleles => 'GT',
        strand  => -1,
        start   => $cds_end - 3,
        end     => $cds_end - 4,
        effects => [qw(frameshift_variant)],
    }, {
        alleles => 'GTAG',
        strand  => -1,
        start   => $cds_end - 3,
        end     => $cds_end - 4,
        effects => [qw(frameshift_variant)],
    }, {
        alleles => '-',
        strand  => -1,
        start   => $cds_end - 3,
        end     => $cds_end - 3,
        effects => [qw(frameshift_variant)],
    }, {
        alleles => '-',
        strand  => -1,
        start   => $cds_end - 4,
        end     => $cds_end - 3,
        effects => [qw(frameshift_variant)],
    }, {
        alleles => '-',
        strand  => -1,
        start   => $cds_end - 6,
        end     => $cds_end - 3,
        effects => [qw(frameshift_variant)],
    }, {
        alleles => 'G',
        strand  => -1,
        start   => $cds_start + 2,
        end     => $cds_start + 2,
        effects => [qw(stop_lost)],
    }, {
        alleles => 'A',
        strand  => -1,
        start   => $cds_start + 1,
        end     => $cds_start + 1,
        effects => [qw(stop_retained_variant)],
    }, {
        alleles => 'C',
        strand  => -1,
        start   => $cds_start,
        end     => $cds_start,
        effects => [qw(stop_lost)],
    }, {
        alleles => 'AAG',
        strand  => -1,
        start   => $cds_start + 2,
        end     => $cds_start + 1,
        effects => [qw(inframe_insertion stop_retained_variant)],
    }, {
        alleles => '-',
        strand  => -1,
        start   => $cds_start,
        end     => $cds_start + 2,
        no_shift => 0,
        effects => [qw(3_prime_UTR_variant coding_sequence_variant)], 
        ## changed for shifting code. Different result is given here than in regular VEP because the transcript
        ## used for the tests is no longer in the gene set, and has the cds_end_NF attribute attached, preventing
        ## overlap_stop_codon from correctly flagging. Test will be updated.
    }, {
        alleles => 'TAA',
        strand  => -1,
        start   => $cds_start,
        end     => $cds_start + 2,
        effects => [qw(stop_retained_variant)],
    }, {
        alleles => 'GGG',
        strand  => -1,
        start   => $cds_start,
        end     => $cds_start + 2,
        effects => [qw(stop_lost)],
    },
    # ============================================================================
    # Bug 2 Reverse Strand Tests: Frameshift indels at stop codon (Issue #1710)
    # Stop codon on reverse strand: $cds_start to $cds_start+2
    # Insertions use start > end (e.g., start=$cds_start+1, end=$cds_start)
    # ============================================================================
    {
        comment => 'Bug2 Reverse: 1bp insertion at stop = frameshift+stop_lost',
        alleles => 'A',
        strand  => -1,
        start   => $cds_start + 1,
        end     => $cds_start,
        effects => [qw(frameshift_variant stop_lost)],
    }, {
        comment => 'Bug2 Reverse: 1bp insertion between 2nd-3rd base = frameshift+stop_lost',
        alleles => 'A',
        strand  => -1,
        start   => $cds_start + 2,
        end     => $cds_start + 1,
        effects => [qw(frameshift_variant stop_lost)],
    }, {
        comment => 'Bug2 Reverse: 2bp insertion at stop = frameshift+stop_lost',
        alleles => 'AT',
        strand  => -1,
        start   => $cds_start + 2,
        end     => $cds_start + 1,
        effects => [qw(frameshift_variant stop_lost)],
    }, {
        comment => 'Bug2 Reverse: 4bp insertion at stop = frameshift+stop_lost',
        alleles => 'ATCG',
        strand  => -1,
        start   => $cds_start + 2,
        end     => $cds_start + 1,
        effects => [qw(frameshift_variant stop_lost)],
    }, {
        comment => 'Bug2 Reverse: 1bp deletion from stop = frameshift+stop_lost',
        alleles => '-',
        strand  => -1,
        start   => $cds_start + 1,
        end     => $cds_start + 1,
        effects => [qw(frameshift_variant stop_lost)],
    }, {
        comment => 'Bug2 Reverse: 2bp deletion from stop = frameshift+stop_lost',
        alleles => '-',
        strand  => -1,
        start   => $cds_start + 1,
        end     => $cds_start + 2,
        effects => [qw(frameshift_variant stop_lost)],
    }, {
        # 3bp insertion at stop codon on reverse strand
        # TGA (stop) + ATG insertion = TATGGA which codes for YG (Tyrosine, Glycine)
        # The stop codon is LOST, not retained - no * in alternate peptide
        # VEP behavior: Only stop_lost is returned (not also inframe_insertion).
        # When a stop codon is lost, VEP doesn't additionally report inframe_insertion.
        comment => 'Bug2 Reverse: 3bp insertion at stop = stop_lost (stop codon disrupted)',
        alleles => 'ATG',
        strand  => -1,
        start   => $cds_start + 2,
        end     => $cds_start + 1,
        effects => [qw(stop_lost)],
    }, {
        comment => 'a wierd allele string',
        alleles => 'HGMD_MUTATION',
        start   => $cds_start + 10,
        end     => $cds_start + 11,
        effects => [qw(coding_sequence_variant)],
    }, {
        comment => 'an ambiguous allele string',
        alleles => 'W',
        start   => $cds_start + 10,
        end     => $cds_start + 10,
        effects => [qw(coding_sequence_variant)],
    }, {
        comment => 'an ambiguous insertion',
        alleles => 'W',
        start   => $cds_start + 10,
        end     => $cds_start + 9,
        effects => [qw(frameshift_variant)],
    }, {
        comment => 'a specified length insertion',
        alleles => '2 BP INSERTION',
        start   => $cds_start + 10,
        end     => $cds_start + 9,
        effects => [qw(frameshift_variant)],
    }, {
        comment => 'an inframe specified length insertion',
        alleles => '3 BP INSERTION',
        start   => $cds_start + 10,
        end     => $cds_start + 9,
        effects => [qw(coding_sequence_variant)],
    }, {
        comment => 'delete the last codon of an exon',
        alleles => '-',
        start   => $intron_end + 1,
        end     => $intron_end + 3,
        effects => [qw(inframe_deletion splice_region_variant)],
    }, 


    # check the complex calls
    
    {
        alleles => '-',
        start   => $intron_end - 2,
        end     => $intron_end + 3,
        effects => [qw(splice_donor_variant splice_donor_region_variant coding_sequence_variant intron_variant)],
    }, {
        alleles => '-',
        start   => $intron_start - 3,
        end     => $intron_start + 2,
        effects => [qw( splice_acceptor_variant coding_sequence_variant intron_variant)],
    }, {
        alleles => '-',
        start   => $cds_end - 2,
        end     => $cds_end + 3,
        effects => [qw( 5_prime_UTR_variant coding_sequence_variant)],
    },  {
        alleles => '-',
        start   => $cds_start - 3,
        end     => $cds_start + 2,
        effects => [qw( 3_prime_UTR_variant coding_sequence_variant)],
    },  


];

# a forward strand transcript with an intron in the UTR

my $t3 = $ta->fetch_by_stable_id('ENST00000530893');

$transcript_tests->{$t3->stable_id}->{transcript} = $t3;

$first_intron = $t3->get_all_Introns->[0];

$intron_start = $first_intron->seq_region_start;
$intron_end   = $first_intron->seq_region_end;

$cds_start = $t3->coding_region_start;

$transcript_tests->{$t3->stable_id}->{tests} = [
    {
        start   => $intron_start - 1,
        end     => $intron_start - 1,
        effects => [qw(splice_region_variant 5_prime_UTR_variant)],
    }, {
        start   => $intron_start,
        end     => $intron_start,
        effects => [qw(splice_donor_variant)],
    }, {
        start   => $intron_start+1,
        end     => $intron_start+1,
        effects => [qw(splice_donor_variant)],
    }, {
        start   => $intron_start+2,
        end     => $intron_start+2,
        effects => [qw(splice_donor_region_variant intron_variant)],
    }, {
        start   => $intron_end - 2,
        end     => $intron_end - 2,
        effects => [qw(splice_region_variant splice_polypyrimidine_tract_variant intron_variant)],
    }, {
        start   => $intron_end - 1,
        end     => $intron_end - 1,
        effects => [qw(splice_acceptor_variant)],
    }, {
        start   => $intron_end,
        end     => $intron_end,
        effects => [qw(splice_acceptor_variant)],
    }, {
        start   => $intron_end + 1,
        end     => $intron_end + 1,
        effects => [qw(splice_region_variant 5_prime_UTR_variant)],
    }, {
        comment => 'a variation with an incorrect reference allele',
        ref     => 'G',
        alleles => 'T',
        start   => $cds_start + 3,
        end     => $cds_start + 3,
        effects => [qw(missense_variant missense_variant)],
    }, {
        comment => 'a TA(3) allele string',
        alleles => 'TA(3)',
        start   => $intron_start + 10,
        end     => $intron_start + 15,
        effects => [qw(intron_variant)],
    },
];

# a forward strand NMD transcript with an intron in the 3 prime UTR

my $nmd_t = $ta->fetch_by_stable_id('ENST00000470094');

$transcript_tests->{$nmd_t->stable_id}->{transcript} = $nmd_t;

my @introns = @{ $nmd_t->get_all_Introns };

my $last_intron = pop @introns;

$intron_start = $last_intron->seq_region_start;
$intron_end   = $last_intron->seq_region_end;

$transcript_tests->{$nmd_t->stable_id}->{tests} = [
    {
        start   => $intron_start - 1,
        end     => $intron_start - 1,
        effects => [qw(splice_region_variant 3_prime_UTR_variant NMD_transcript_variant)],
    }, {
        start   => $intron_start + 1,
        end     => $intron_start + 1,
        effects => [qw(splice_donor_variant NMD_transcript_variant)],
    }, {
        start   => $intron_end + 1,
        end     => $intron_end + 1,
        effects => [qw(splice_region_variant 3_prime_UTR_variant NMD_transcript_variant)],
    }, 
];

# a miRNA transcript

my $mirna = $ta->fetch_by_stable_id('ENST00000408781');

$transcript_tests->{$mirna->stable_id}->{transcript} = $mirna;

$t_start = $mirna->seq_region_start;
$t_end   = $mirna->seq_region_end;

$transcript_tests->{$mirna->stable_id}->{tests} = [
    {
        start   => $t_start,
        end     => $t_start,
        effects => [qw(non_coding_transcript_exon_variant)],
    }, {
        start   => $t_start + 40,
        end     => $t_start + 40,
        effects => [qw(mature_miRNA_variant)],
    }, 
];

# a forward strand transcript with a partial stop codon

my $t4 = $ta->fetch_by_stable_id('ENST00000450073');

$transcript_tests->{$t4->stable_id}->{transcript} = $t4;

$cds_start = $t4->coding_region_start;
$cds_end   = $t4->coding_region_end;

$transcript_tests->{$t4->stable_id}->{tests} = [
    {
        start   => $cds_end,
        end     => $cds_end,
        effects => [qw(incomplete_terminal_codon_variant coding_sequence_variant)],
    }, 
    {
        start   => $cds_end,
        end     => $cds_end,
        alleles => '-',
        effects => [qw(incomplete_terminal_codon_variant coding_sequence_variant)],
    }, 
];

# transcripts with frameshift introns

# a non-coding pseudogene with a frameshift intron
my $nc_fs_t = $ta->fetch_by_stable_id('ENST00000438775');

$transcript_tests->{$nc_fs_t->stable_id}->{transcript} = $nc_fs_t;

my $fs_intron_start = $nc_fs_t->get_all_Exons->[0]->end+1;

$transcript_tests->{$nc_fs_t->stable_id}->{tests} = [
    {
        comment => "a non-coding transcript with a frameshift intron",
        start   => $fs_intron_start,
        end     => $fs_intron_start,
        effects => [qw(non_coding_transcript_variant)],
    }, 
];

# a coding transcript with a frameshift intron in the CDS
my $c_fs_t = $ta->fetch_by_stable_id('ENST00000392535');

$transcript_tests->{$c_fs_t->stable_id}->{transcript} = $c_fs_t;

$fs_intron_start = $c_fs_t->get_all_Exons->[0]->end+1;

$transcript_tests->{$c_fs_t->stable_id}->{tests} = [
    {
        comment => "a transcript with a frameshift intron in the CDS",
        start   => $fs_intron_start,
        end     => $fs_intron_start,
        effects => [qw(coding_sequence_variant)],
    }, 
];

# a transcript with a selenocysteine edit
my $sc_se_t = $ta->fetch_by_stable_id('ENST00000380903');

$transcript_tests->{$sc_se_t->stable_id}->{transcript} = $sc_se_t;

$transcript_tests->{$sc_se_t->stable_id}->{tests} = [
    {
        comment => "a transcript with a selenocysteine seqEdit",
        alleles => 'C',
        start   => 50655788,
        end     => 50655788,
        effects => [qw(missense_variant)],
    }, 
];

# a transcript with a misc amino acid edit altered to start
my $aa_se_t = $ta->fetch_by_stable_id('ENST00000295641');

$transcript_tests->{$aa_se_t->stable_id}->{transcript} = $aa_se_t;

$transcript_tests->{$aa_se_t->stable_id}->{tests} = [
    {
        comment => "a transcript with a misc amino acid seqEdit",
        alleles => 'T',
        start   => 220462640,
        end     => 220462640,
        effects => [qw(start_retained_variant)],
    }, 
];

# a transcript with incomplete 5' CDS
my $incomplete_cds_t = $ta->fetch_by_stable_id('ENST00000452863');
$transcript_tests->{$incomplete_cds_t->stable_id}->{transcript} = $incomplete_cds_t;

$transcript_tests->{$incomplete_cds_t->stable_id}->{tests} = [
    {
        comment => "a transcript with incomplete 5' CDS",
        alleles => 'T',
        start   =>  32456435,
        end     =>  32456435,
        effects => [qw(missense_variant)],
    }, 
];

my $transcript = $transcript_tests->{$incomplete_cds_t->stable_id}->{transcript};
my $vf = Bio::EnsEMBL::Variation::VariationFeature->new(
    -start          => 220462640,
    -end            => 220462640,
    -strand         => 1,
    -slice          => $transcript->slice,
    -allele_string  => 'A/G',
    -variation_name => 'test_start_retained',
);

my $tv = Bio::EnsEMBL::Variation::TranscriptVariation->new(
    -variation_feature  => $vf,
    -transcript         => $transcript,
);

my $tva = $tv->get_all_alternate_BaseVariationFeatureOverlapAlleles();

my $start_retained = Bio::EnsEMBL::Variation::Utils::VariationEffect::start_retained_variant($tva->[0]);
is($start_retained, 0, 'start_retained works with no $bvfo & $bvf');

my $stop_retained = Bio::EnsEMBL::Variation::Utils::VariationEffect::stop_retained($tva->[0]);
is($stop_retained, undef, 'stop_retained works with no $bvfo & $bvf');

my $coding_unknown = Bio::EnsEMBL::Variation::Utils::VariationEffect::coding_unknown($tva->[0]);
is($coding_unknown, 0, 'coding_unknown works with no $bvfo & $bvf');

my $bvfo = $tva->[0]->base_variation_feature_overlap;
my $bvf = $bvfo->base_variation_feature;
$bvf->{allele_string} = 'COSMIC_MUTATION';
my $start_retained_cosmic = Bio::EnsEMBL::Variation::Utils::VariationEffect::start_retained_variant($tva->[0], 0, $bvfo, $bvf);
is($start_retained_cosmic, 0, 'start_retained retuns 0 with COSMIC');

delete($tva->[0]->{_predicate_cache}->{stop_retained});
my $stop_retained_cosmic = Bio::EnsEMBL::Variation::Utils::VariationEffect::stop_retained($tva->[0], 0, $bvfo, $bvf);
is($stop_retained_cosmic, 0, 'stop_retained returns 0 with COSMIC');

my $coding_unknown_cosmic = Bio::EnsEMBL::Variation::Utils::VariationEffect::coding_unknown($tva->[0], 0, $bvfo, $bvf);
is($coding_unknown_cosmic, 0, 'coding_unknown returns 0 with COSMIC');

my $vf_cosmic = Bio::EnsEMBL::Variation::VariationFeature->new(
    -start          => 20462640,
    -end            => 20462640,
    -strand         => 1,
    -slice          => $transcript->slice,
    -allele_string  => '-/COSMIC_MUTATION',
    -variation_name => 'test_cosmic_shift',
);

$vf_cosmic->{class_display_term} = 'insertion';

my $tv_cosmic = Bio::EnsEMBL::Variation::TranscriptVariation->new(
    -variation_feature  => $vf_cosmic,
    -transcript         => $transcript,
);

## Check that COSMIC_MUTATIONS are not shifted
my $tva_cosmic = $tv_cosmic->get_all_alternate_BaseVariationFeatureOverlapAlleles();
my $bvfo_cosmic = $tva_cosmic->[0]->base_variation_feature_overlap;
my $bvf_cosmic = $bvfo_cosmic->base_variation_feature;
$bvf_cosmic->{tva_shift_hashes} = [];
$tva_cosmic->[0]->_return_3prime;
is($tva_cosmic->[0]->{shift_hash}, undef, 'COSMIC_MUTATIONs has no shift hash');

$bvf_cosmic->{allele_string} = '-/G';
$tva_cosmic->[0]->{allele_string} = '-/G';
$tva_cosmic->[0]->_return_3prime;
is(defined($tva_cosmic->[0]->{shift_hash}), 1, 'non-COSMIC_MUTATIONs has shift hash');

my $test_count = 1;

my $def_strand  = 1;

my $reverse = 0;

for my $stable_id (keys %$transcript_tests) {
    
    my $tran = $transcript_tests->{$stable_id}->{transcript};

    for my $test (@{ $transcript_tests->{$stable_id}->{tests} }) {

        my $ref = $test->{ref} || $tran->slice->subseq($test->{start}, $test->{end}, $test->{strand});
       
        $ref = '-' unless $ref;

        unless ($test->{alleles}) {
            my $alt = $ref;
            reverse_comp(\$alt);
            $test->{alleles} = $alt;
        }

        $DB::single = 1 if $test->{debug};

        $test->{strand} = $def_strand unless defined $test->{strand};

        my $no_shift = $test->{no_shift};

        my $allele_string = $ref.'/'.$test->{alleles};

        my $vf = Bio::EnsEMBL::Variation::VariationFeature->new(
            -start          => $test->{start},
            -end            => $test->{end},
            -strand         => $test->{strand},
            -slice          => $tran->slice,
            -allele_string  => $ref.'/'.$test->{alleles},
            -variation_name => 'test'.$test_count,
        );

        my $tv = Bio::EnsEMBL::Variation::TranscriptVariation->new(
            -variation_feature  => $vf,
            -transcript         => $tran,
            -no_shift		=> $no_shift,
	);

        warn "# alleles: $allele_string\n";
        warn '# codons: ', $tv->codons, "\n" if $tv->codons;
        warn '# peptides: ', $tv->pep_allele_string, "\n" if $tv->pep_allele_string;

        my @effects = map {
            map { $_->SO_term } @{ $_->get_all_OverlapConsequences }
        } @{ $tv->get_all_alternate_TranscriptVariationAlleles };

        my $comment = $test->{comment} || (join ',', @{ $test->{effects} }) || 'no effect';


        my $strand = $tv->transcript->strand;

        # sort so that the order doesn't matter
        is_deeply( [sort @effects], [sort @{ $test->{effects} }], "VF $test_count (strand $strand): $comment") 
            || (diag "Actually got: ", explain \@effects, "but expected: ", explain \@{ $test->{effects}} ) ||  
            print join(",",(sort @{ $test->{effects}} )) ."\tgot\t".  join(",",(sort @effects)) . "\t";

        if($DEBUG ==1){
          print $tv->hgvs_genomic()->{$test->{alleles}} ."\t" ;
          print $tv->hgvs_transcript()->{$test->{alleles}} ."\t" if defined  $tv->hgvs_transcript()->{$test->{alleles}};
          print $tv->hgvs_protein()->{$test->{alleles}} if defined  $tv->hgvs_protein()->{$test->{alleles}};
          print "\n" ;
        }

        if (my $expected_pep_alleles = $test->{pep_alleles}) {
            is(
                $tv->pep_allele_string, 
                $expected_pep_alleles, 
                "peptide allele string is correct (expected $expected_pep_alleles)"
            ) || die;
        }

        $test_count++;
    }
}


done_testing();

