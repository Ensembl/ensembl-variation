package Bio::EnsEMBL::Variation::Utils::VariationEffect;

use strict;
use warnings;

use Data::Dumper;

#use Inline C => <<'END_C';
#
#int overlap (int f1_start, int f1_end, int f2_start, int f2_end) {
#    return (f1_end >= f2_start && f1_start <= f2_end);
#}
#
#END_C

sub overlap {
    my ( $f1_start, $f1_end, $f2_start, $f2_end ) = @_;
    return ($f1_end >= $f2_start and $f1_start <= $f2_end);
}

sub overlaps_transcript {
    my ($vf, $tran) = @_;
    
    return 0 unless $tran->isa('Bio::EnsEMBL::Transcript');
    
    return overlap($vf->start, $vf->end, $tran->start - 5000, $tran->end + 5000);
}

sub within_transcript {
    my $tva     = shift;
    my $vf      = $tva->variation_feature;
    my $tran    = $tva->transcript;
    
    return overlap($vf->start, $vf->end, $tran->start, $tran->end);
}

sub _before_start {
    my ($vf, $tran, $dist) = @_;
    return ( ($vf->end >= ($tran->start - $dist)) and ($vf->end < $tran->start) );
}

sub _after_end {
    my ($vf, $tran, $dist) = @_;
    return ( ($vf->start <= ($tran->end + $dist)) and ($vf->start > $tran->end) );
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
    my $tva    = shift;
    my $vf      = $tva->variation_feature;
    my $tran    = $tva->transcript;

    return _upstream($vf, $tran, 5000);
}

sub downstream_5KB {
    my $tva     = shift;
    my $vf      = $tva->variation_feature;
    my $tran    = $tva->transcript; 

    return _downstream($vf, $tran, 5000);
}

sub upstream_2KB {
    my $tva     = shift;
    my $vf      = $tva->variation_feature;
    my $tran    = $tva->transcript; 

    return _upstream($vf, $tran, 2000);
}

sub downstream_2KB {
    my $tva     = shift;
    my $vf      = $tva->variation_feature;
    my $tran    = $tva->transcript; 

    return _downstream($vf, $tran, 2000);
}

sub downstream_500B {
    my $tva     = shift;
    my $vf      = $tva->variation_feature;
    my $tran    = $tva->transcript; 

    return _downstream($vf, $tran, 500);
}

sub within_nmd_transcript {
    my $tva     = shift;
    my $tran    = $tva->transcript; 

    return ( within_transcript($tva) and ($tran->biotype eq 'nonsense_mediated_decay') );
}

sub within_non_coding_gene {
    my $tva     = shift;
    my $tran    = $tva->transcript;
    
    return ( within_transcript($tva) and (not $tran->translation) );
}

sub within_miRNA {
    my $tva     = shift;
    my $tv      = $tva->transcript_variation;
    my $vf      = $tva->variation_feature;
    my $tran    = $tva->transcript;
    
    if ($tran->biotype eq 'miRNA') {
        my ($attribute) = @{ $tran->get_all_Attributes('miRNA') };
        
        if (defined $attribute && $attribute->value =~ /(\d+)-(\d+)/) { 
            for my $coord ($tv->mapper->cdna2genomic($1, $2, $tran->strand)) {
                if ($coord->isa('Bio::EnsEMBL::Mapper::Coordinate')) {
                    if (overlap($vf->start, $vf->end, $coord->start, $coord->end)) {
                        return 1;
                    }
                }
            }
        }
    }
    
    return 0;
}

sub donor_splice_site {
    my $tva     = shift;
    my $tran    = $tva->transcript;
    
    return $tran->strand == 1 ? 
        $tva->transcript_variation->intron_effects->{start_splice_site} :
        $tva->transcript_variation->intron_effects->{end_splice_site};
}

sub acceptor_splice_site {
    my $tva     = shift;
    my $tran    = $tva->transcript;
    
    return $tran->strand == 1 ? 
        $tva->transcript_variation->intron_effects->{end_splice_site} :
        $tva->transcript_variation->intron_effects->{start_splice_site};
}

sub essential_splice_site {
    my $tva = shift;
    
    return ( acceptor_splice_site($tva) or donor_splice_site($tva) );
}

sub splice_region {
    my $tva    = shift;

    return $tva->transcript_variation->intron_effects->{splice_region};
}

sub within_intron {
    my $tva    = shift;

    return $tva->transcript_variation->intron_effects->{intronic};
}

sub within_coding_region {
    my $tva     = shift;
    my $vf      = $tva->variation_feature;
    my $tran    = $tva->transcript;
    
    return 0 unless $tran->translation;
    
    return overlap($vf->start, $vf->end, $tran->coding_region_start, $tran->coding_region_end);
}

sub _before_coding {
    my ($vf, $tran) = @_;
    return 0 unless defined $tran->translation;
    return overlap($vf->start, $vf->end, $tran->start, $tran->coding_region_start-1);    
}

sub _after_coding {
    my ($vf, $tran) = @_;
    return 0 unless defined $tran->translation;
    return overlap($vf->start, $vf->end, $tran->coding_region_end+1, $tran->end);    
}

sub within_5_prime_utr {
    my $tva     = shift;
    my $vf      = $tva->variation_feature;
    my $tran    = $tva->transcript;
    
    my $five_prime_of_coding = 
        $tran->strand == 1 ? 
            _before_coding($vf, $tran) : 
            _after_coding($vf, $tran);
    
    return ( $five_prime_of_coding and (not within_intron($tva)) );
}

sub within_3_prime_utr {
    my $tva     = shift;
    my $vf      = $tva->variation_feature;
    my $tran    = $tva->transcript; 
    
    my $three_prime_of_coding = 
        $tran->strand == 1 ? 
            _after_coding($vf, $tran) : 
            _before_coding($vf, $tran);
    
    return ( $three_prime_of_coding and (not within_intron($tva)) );
}

sub complex_indel {
    my $tva     = shift;
    my $vf      = $tva->variation_feature;
    
    return 0 unless $vf->var_class eq 'in-del';

    return @{ $tva->transcript_variation->cds_coords } > 1;
}

sub synonymous_coding {
    my $tva = shift;
    my $tv  = $tva->transcript_variation;

    return 0 unless $tva->affects_cds;
    
    my $alt_seq = $tva->feature_seq;
    my $ref_seq = $tv->reference_allele->feature_seq;
    
    return 0 if ($alt_seq eq '-' or $ref_seq eq '-');
    
    my $alt_aa = $tva->amino_acid;
    my $ref_aa = $tv->reference_allele->amino_acid;

    return ( $alt_aa eq $ref_aa );
}

sub non_synonymous_coding {
    my $tva = shift;
    my $tv  = $tva->transcript_variation;

    return 0 unless $tva->affects_cds;
    
    my $alt_seq = $tva->feature_seq;
    my $ref_seq = $tv->reference_allele->feature_seq;
    
    return 0 if ($alt_seq eq '-' or $ref_seq eq '-');

    my $alt_aa = $tva->amino_acid;
    my $ref_aa = $tv->reference_allele->amino_acid;
    
    return ( $alt_aa ne $ref_aa );
}

sub stop_gained {
    my $tva = shift;
    my $tv  = $tva->transcript_variation;

    return 0 unless $tva->affects_cds;
    
    my $alt_aa = $tva->amino_acid;
    my $ref_aa = $tv->reference_allele->amino_acid;

    return ( ($alt_aa =~ /\*/) and ($ref_aa !~ /\*/) );
}

sub stop_lost {
    my $tva = shift;
    my $tv  = $tva->transcript_variation;

    return 0 unless $tva->affects_cds;
    
    my $alt_aa = $tva->amino_acid;
    my $ref_aa = $tv->reference_allele->amino_acid;

    return ( ($alt_aa !~ /\*/) and ($ref_aa =~ /\*/) );
}

sub frameshift {
    my $tva = shift;

    return 0 unless $tva->affects_cds;
    
    return 0 if partial_codon($tva);

    my $tv = $tva->transcript_variation;
    my $var_len = $tv->cds_end - $tv->cds_start + 1;

    my $seq = $tva->feature_seq;
    
    $seq = '' if $seq eq '-';

    return abs( length($seq) - $var_len ) % 3;
}

sub partial_codon {
    my $tva = shift;
    
    return 0 unless $tva->affects_cds;

    my $tv = $tva->transcript_variation;

    my $cds_length = length $tv->translateable_seq;

    my $codon_cds_start = ($tv->pep_start * 3) - 2;

    my $last_codon_length = $cds_length - ($codon_cds_start - 1);
    
    return ( $last_codon_length < 3 and $last_codon_length > 0 );
}

sub within_coding_frameshift_intron {
    my $tva = shift;
    
    return (within_coding_region($tva) and 
        $tva->transcript_variation->intron_effects->{within_frameshift_intron});
}

1;

