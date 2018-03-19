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
        alleles => 'A',
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
        effects => [qw(splice_region_variant intron_variant)],
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
        effects => [qw(intron_variant)],
    }, {
        start   => $intron_end - 7,
        end     => $intron_end - 7,
        effects => [qw(splice_region_variant intron_variant)],
    }, {
        start   => $intron_end - 2,
        end     => $intron_end - 2,
        effects => [qw(splice_region_variant intron_variant)],
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
        alleles => 'A',
        start   => $intron_end-1,
        end     => $intron_end-2,
        effects => [qw(splice_region_variant intron_variant)],
    }, {
        comment => 'insertion between last bases 7 & 8 is still splice_region',
        alleles => 'A',
        start   => $intron_end-6,
        end     => $intron_end-7,
        effects => [qw(splice_region_variant intron_variant)],
    }, {
        comment => 'insertion between last bases 8 & 9 is just an intron_variant',
        alleles => 'A',
        start   => $intron_end-7,
        end     => $intron_end-8,
        effects => [qw(intron_variant)],
    }, {
        comment => 'whole exon deletion',
        alleles => '-',
        start   => $exon_start-10,
        end     => $exon_end+10,
        effects => [qw(intron_variant splice_acceptor_variant splice_donor_variant coding_sequence_variant )],
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
        effects => [qw(inframe_insertion)],
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
        effects => [qw(5_prime_UTR_variant start_lost)],
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
    }, {
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
        comment => 'delete the last codon of an exon',
        alleles => '-',
        start   => $intron_start-3,
        end     => $intron_start-1,
        effects => [qw(inframe_deletion splice_region_variant)],
    }, 
    

    # check the complex calls
   
    {
        alleles => '-',
        start   => $intron_start-3,
        end     => $intron_start+2,
        effects => [qw( splice_donor_variant coding_sequence_variant intron_variant)],
    }, {
        alleles => '-',
        start   => $intron_end-2,
        end     => $intron_end+3,
        effects => [qw( splice_acceptor_variant coding_sequence_variant  intron_variant)],
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
        comment => 'deletion overlapping STOP and 3\' UTR, stop retained',
        alleles => '-',
        start   => $cds_end-1,
        end     => $cds_end+1,
        effects => [qw( 3_prime_UTR_variant stop_retained_variant)],
    }, {
        comment => 'deletion overlapping STOP and 3\' UTR, stop retained, different codon',
        alleles => '-',
        start   => $cds_end-1,
        end     => $cds_end+2,
        effects => [qw( 3_prime_UTR_variant stop_retained_variant)],
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
        alleles => 'A',
        start   => $t_end + 1,
        end     => $t_end,
        effects => [ qw(upstream_gene_variant) ],
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
        effects => [qw(splice_region_variant intron_variant)],
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
        alleles => 'A',
        start   => $intron_end,
        end     => $intron_end - 1,
        effects => [qw(splice_donor_variant)],
    }, {
        comment => 'insertion between bases 2 & 3 of an intron is splice_region',
        alleles => 'A',
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
        effects => [qw(intron_variant)],
    }, {
        start   => $intron_start + 7,
        end     => $intron_start + 7,
        effects => [qw(splice_region_variant intron_variant)],
    }, {
        start   => $intron_start + 2,
        end     => $intron_start + 2,
        effects => [qw(splice_region_variant intron_variant)],
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
        effects => [qw(splice_region_variant intron_variant)],
    }, {
        comment => 'insertion between last bases 7 & 8 is still splice_region',
        alleles => 'A',
        start   => $intron_start + 7,
        end     => $intron_start + 6,
        effects => [qw(splice_region_variant intron_variant)],
    }, {
        comment => 'insertion between last bases 8 & 9 is just an intron_variant',
        alleles => 'A',
        start   => $intron_start + 8,
        end     => $intron_start + 7,
        effects => [qw(intron_variant)],
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
        effects => [qw(stop_lost inframe_deletion)],
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
        effects => [qw(splice_donor_variant coding_sequence_variant intron_variant)],
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
        effects => [qw(splice_region_variant intron_variant)],
    }, {
        start   => $intron_end - 2,
        end     => $intron_end - 2,
        effects => [qw(splice_region_variant intron_variant)],
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

# a transcript with a misc amino acid edit
my $aa_se_t = $ta->fetch_by_stable_id('ENST00000295641');

$transcript_tests->{$aa_se_t->stable_id}->{transcript} = $aa_se_t;

$transcript_tests->{$aa_se_t->stable_id}->{tests} = [
    {
        comment => "a transcript with a misc amino acid seqEdit",
        alleles => 'T',
        start   => 220462640,
        end     => 220462640,
        effects => [qw(synonymous_variant)],
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

