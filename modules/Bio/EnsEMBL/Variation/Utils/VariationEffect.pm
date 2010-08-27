package Bio::EnsEMBL::Variation::Utils::VariationEffect;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::VariationFeatureOverlap;
use Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele;
use Bio::EnsEMBL::Variation::OverlapConsequence;

use vars qw(@ISA @EXPORT_OK);

@ISA = qw(Exporter);

@EXPORT_OK = qw(&transcript_effect &regulatory_region_effect &get_all_variation_effects);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);

use Data::Dumper;

sub overlap {
    my ( $f1_start, $f1_end, $f2_start, $f2_end ) = @_;
    return ($f1_end >= $f2_start and $f1_start <= $f2_end);
}

sub overlaps_transcript {
    my ($vf, $tran) = @_;
    return overlap($vf->start, $vf->end, $tran->start - 5000, $tran->end + 5000);
}

sub within_transcript {
    my $vfoa    = shift;
    my $vf      = $vfoa->variation_feature_overlap->variation_feature;
    my $tran    = $vfoa->variation_feature_overlap->feature; 
    return overlap($vf->start, $vf->end, $tran->start, $tran->end);
}

sub _before_start {
    my ($vf, $tran, $dist) = @_;
    return ($vf->end >= ($tran->start - $dist) and $vf->end < $tran->start);
}

sub _after_end {
    my ($vf, $tran, $dist) = @_;
    return ($vf->start <= ($tran->end + $dist) and $vf->start > $tran->end);
}

sub _upstream {
    my ($vf, $tran, $dist) = @_;
    return $tran->strand == 1 ? _before_start($vf, $tran, $dist) : _after_end($vf, $tran, $dist);
}

sub _downstream {
    my ($vf, $tran, $allele, $dist) = @_;
    return $tran->strand == 1 ? _after_end($vf, $tran, $allele, $dist) : _before_start($vf, $tran, $allele, $dist);
}

sub upstream_5KB {
    my $vfoa    = shift;
    my $vf      = $vfoa->variation_feature_overlap->variation_feature;
    my $tran    = $vfoa->variation_feature_overlap->feature; 
   
    return 0 unless $tran->isa('Bio::EnsEMBL::Transcript');

    return _upstream($vf, $tran, 5000);
}

sub downstream_5KB {
    my $vfoa    = shift;
    my $vf      = $vfoa->variation_feature_overlap->variation_feature;
    my $tran    = $vfoa->variation_feature_overlap->feature; 

    return _downstream($vf, $tran, 5000);
}

sub upstream_2KB {
    my $vfoa    = shift;
    my $vf      = $vfoa->variation_feature_overlap->variation_feature;
    my $tran    = $vfoa->variation_feature_overlap->feature; 

    return _upstream($vf, $tran, 2000);
}

sub downstream_2KB {
    my $vfoa    = shift;
    my $vf      = $vfoa->variation_feature_overlap->variation_feature;
    my $tran    = $vfoa->variation_feature_overlap->feature; 

    return _downstream($vf, $tran, 2000);
}

sub within_nmd_transcript {
    my $vfoa    = shift;
    my $vf      = $vfoa->variation_feature_overlap->variation_feature;
    my $tran    = $vfoa->variation_feature_overlap->feature; 

    return (within_transcript($vfoa) and $tran->biotype eq 'nonsense_mediated_decay');
}

sub within_non_coding_gene {
    my $vfoa    = shift;
    my $vf      = $vfoa->variation_feature_overlap->variation_feature;
    my $tran    = $vfoa->variation_feature_overlap->feature;
    
    return 0 if within_miRNA($vfoa);
    
    return (within_transcript($vfoa) and (not $tran->translation));
}

sub within_miRNA {
    my $vfoa    = shift;
    my $tran    = $vfoa->variation_feature_overlap->feature;
    
    return 0;
}

sub _start_splice_site {
    my ($vf, $tran) = @_;

    for my $intron (@{ $tran->get_all_Introns }) {
        if (overlap($vf->start, $vf->end, $intron->start, $intron->start+2)) {
            return 1;
        }
    }

    return 0;
}

sub _end_splice_site {
    my ($vf, $tran) = @_;
    
    for my $intron (@{ $tran->get_all_Introns }) {
        if (overlap($vf->start, $vf->end, $intron->end-2, $intron->end)) {
            return 1;
        }
    }

    return 0;
}

sub donor_splice_site {
    my $vfoa    = shift;
    my $vf      = $vfoa->variation_feature_overlap->variation_feature;
    my $tran    = $vfoa->variation_feature_overlap->feature; 
    
    return $tran->strand == 1 ? _start_splice_site($vf, $tran) : _end_splice_site($vf, $tran);
}

sub acceptor_splice_site {
    my $vfoa    = shift;
    my $vf      = $vfoa->variation_feature_overlap->variation_feature;
    my $tran    = $vfoa->variation_feature_overlap->feature; 
    
    return $tran->strand == 1 ? _end_splice_site($vf, $tran) : _start_splice_site($vf, $tran); 
}

sub essential_splice_site {
    my $vfoa    = shift;
    return (acceptor_splice_site($vfoa) or donor_splice_site($vfoa));
}

sub splice_region {
    my $vfoa    = shift;
    my $vf      = $vfoa->variation_feature_overlap->variation_feature;
    my $tran    = $vfoa->variation_feature_overlap->feature; 

    for my $intron (@{ $tran->get_all_Introns }) {

        if ( overlap($vf->start, $vf->end, $intron->start-3, $intron->start+8) or
             overlap($vf->start, $vf->end, $intron->end-8, $intron->end+3) ) {
            return 1;
        }
    }

    return 0;
}

sub within_intron {
    my $vfoa    = shift;
    my $vf      = $vfoa->variation_feature_overlap->variation_feature;
    my $tran    = $vfoa->variation_feature_overlap->feature; 

    for my $intron (@{ $tran->get_all_Introns }) {
        if (overlap($vf->start, $vf->end, $intron->start, $intron->end)) {
            return 1;
        }
    }

    return 0;
}

sub _start_utr {
    my ($vf, $tran) = @_;
    return 0 unless defined $tran->translation;
    return overlap($vf->start, $vf->end, $tran->start, $tran->coding_region_start-1);    
}

sub _end_utr {
    my ($vf, $tran) = @_;
    return 0 unless defined $tran->translation;
    return overlap($vf->start, $vf->end, $tran->coding_region_end+1, $tran->end);    
}

sub within_5_prime_utr {
    my $vfoa    = shift;
    my $vf      = $vfoa->variation_feature_overlap->variation_feature;
    my $tran    = $vfoa->variation_feature_overlap->feature; 
    
    return $tran->strand == 1 ? _start_utr($vf, $tran) : _end_utr($vf, $tran); 
}

sub within_3_prime_utr {
    my $vfoa    = shift;
    my $vf      = $vfoa->variation_feature_overlap->variation_feature;
    my $tran    = $vfoa->variation_feature_overlap->feature; 
    
    return $tran->strand == 1 ? _end_utr($vf, $tran) : _start_utr($vf, $tran); 
}

sub complex_indel {
    my $vfoa    = shift;
    my $vf      = $vfoa->variation_feature_overlap->variation_feature;
    
    return 0 unless $vf->var_class eq 'in-del';

    return @{ $vfoa->variation_feature_overlap->cds_coords } > 1;
}

sub synonymous_coding {
    my $vfoa = shift;

    return 0 unless $vfoa->affects_cds;
    
    my $alt_aa = $vfoa->aa;
    my $ref_aa = $vfoa->variation_feature_overlap->reference_allele->aa;
    
    print "ref_aa: $ref_aa alt_aa: $alt_aa\n";
    
    return 0 if ($alt_aa eq '-' or $ref_aa eq '-');
    
    return $alt_aa eq $ref_aa;
}

sub non_synonymous_coding {
    my $vfoa = shift;

    return 0 unless $vfoa->affects_cds;

    my $alt_aa = $vfoa->aa;
    my $ref_aa = $vfoa->variation_feature_overlap->reference_allele->aa;
    
    return 0 if ($alt_aa eq '-' or $ref_aa eq '-');
    
    return $alt_aa ne $ref_aa;
}

sub stop_gained {
    my $vfoa = shift;

    return 0 unless $vfoa->affects_cds;
    
    my $alt_aa = $vfoa->aa;
    my $ref_aa = $vfoa->variation_feature_overlap->reference_allele->aa;
    
    return 0 unless $alt_aa and $ref_aa;

    return ($alt_aa =~ /\*/ and $ref_aa !~ /\*/);
}

sub stop_lost {
    my $vfoa = shift;

    return 0 unless $vfoa->affects_cds;
    
    my $alt_aa = $vfoa->aa;
    my $ref_aa = $vfoa->variation_feature_overlap->reference_allele->aa;
    
    return 0 unless $alt_aa and $ref_aa;

    return ($alt_aa !~ /\*/ and $ref_aa =~ /\*/);
}

sub frameshift {
    my $vfoa = shift;

    return 0 unless $vfoa->affects_cds;
    
    return 0 if partial_codon($vfoa);

    my $vfo = $vfoa->variation_feature_overlap;
    my $var_len = $vfo->cds_end - $vfo->cds_start + 1;

    my $seq = $vfoa->seq;
    
    $seq = '' if $seq eq '-';

    return abs( length($seq) - $var_len ) % 3;
}

sub partial_codon {
    my $vfoa = shift;
    
    return 0 unless $vfoa->affects_cds;

    my $vfo = $vfoa->variation_feature_overlap;

    my $cds_length = length($vfo->feature->translateable_seq);

    my $codon_cds_start = ($vfo->pep_start * 3) - 2;

    my $last_codon_length = $cds_length - ($codon_cds_start - 1);
    
    return ($last_codon_length < 3 and $last_codon_length > 0);
}

my $effect_rules = { 
    '5KB_upstream_variant'              => \&upstream_5KB,
    '5KB_downstream_variant'            => \&downstream_5KB,
    '2KB_upstream_variant'              => \&upstream_2KB,
    '2KB_downstream_variant'            => \&downstream_2KB,
    'splice_donor_variant'              => \&donor_splice_site,
    'splice_acceptor_variant'           => \&acceptor_splice_site,
    'splice_site_variant'               => \&essential_splice_site,
    'splice_region_variant'             => \&splice_region,
    'intron_variant'                    => \&within_intron,
    '5_prime_UTR_variant'               => \&within_5_prime_utr,
    '3_prime_UTR_variant'               => \&within_3_prime_utr,
    'complex_change_in_transcript'      => \&complex_indel,
    'synonymous_codon'                  => \&synonymous_coding,
    'non_synonymous_codon'              => \&non_synonymous_coding,
    'stop_gained'                       => \&stop_gained,
    'stop_lost'                         => \&stop_lost,
    'frameshift_variant'                => \&frameshift,
    'incomplete_terminal_codon_variant' => \&partial_codon,
    'NMD_transcript_variant'            => \&within_nmd_transcript,
    'non_coding_RNA_variant'            => \&within_non_coding_gene,
    'mature_miRNA_variant'              => \&within_miRNA,
    
};

my @consequences;

for my $effect (keys %$effect_rules) {
    my $predicate = $effect_rules->{$effect};
    push @consequences, Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
            SO_term     => $effect,
            predicate   => $predicate,
        });
}

sub transcript_effect {

    my ($vf, $tran) = @_;

    my $result = {};

    print "T: ",$tran->start, "-", $tran->end," ",$tran->strand,"\n";
    print "V: ",$vf->start, "-", $vf->end," ",$vf->strand,"\n";
    print $tran->stable_id,"\n";
    print $vf->variation_name,"\n";

    if (overlaps_transcript($vf, $tran)) {

        my $vfo = Bio::EnsEMBL::Variation::VariationFeatureOverlap->new_fast({            
            variation_feature  => $vf,
            feature            => $tran,
            feature_type       => 'transcript',
        });
        
        my $tm = $tran->get_TranscriptMapper;

        $vfo->cdna_coords([$tm->genomic2cdna($vf->start, $vf->end, $tran->strand)]);
        $vfo->cds_coords([$tm->genomic2cds($vf->start, $vf->end, $tran->strand)]);
        $vfo->pep_coords([$tm->genomic2pep($vf->start, $vf->end, $tran->strand)]);

        my @alleles = split /\//, $vf->allele_string;

        # if the strands differ we need to reverse complement the alleles from the allele_string
        unless ($tran->strand == $vf->strand) {
            map { reverse_comp(\$_) } @alleles;
        }

        my $ref_allele = shift @alleles;

        my $ref_vfoa = Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele->new_fast({
            variation_feature_overlap   => $vfo,
            seq                         => $ref_allele,
            is_reference                => 1,
        });

        $vfo->reference_allele($ref_vfoa);
        $vfo->alleles($ref_vfoa);

        for my $allele (@alleles) {
            
            my $vfoa = Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele->new_fast({
                variation_feature_overlap   => $vfo,
                seq                         => $allele,
            });
            $vfoa->calc_consequences(\@consequences);
            $vfo->alleles($vfoa);
        }

        return $vfo;
    }
    else {
        return undef;
    }
}

sub regulatory_region_effect {
    my ($reg_region, $vf) = @_;

}

sub get_all_variation_effects {

    my ($vf) = @_;

    # get all features that overlap this variation
    #
    #

}

1;

