=head1 LICENSE

 Copyright (c) 1999-2011 The European Bioinformatics Institute and
 Genome Research Limited.  All rights reserved.

 This software is distributed under a modified Apache license.
 For license details, please see

   http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <dev@ensembl.org>.

 Questions may also be sent to the Ensembl help desk at
 <helpdesk@ensembl.org>.

=cut

package Bio::EnsEMBL::Variation::Utils::VariationEffect;

use strict;
use warnings;

use base qw(Exporter);

our @EXPORT_OK = qw(overlap affects_peptide);

sub overlap {
    my ( $f1_start, $f1_end, $f2_start, $f2_end ) = @_;
    return ( ($f1_end >= $f2_start) and ($f1_start <= $f2_end) );
}

sub overlaps_transcript {
    my ($vf, $tran) = @_;
    
    return 0 unless $tran->isa('Bio::EnsEMBL::Transcript');
    
    return overlap($vf->seq_region_start, $vf->seq_region_end,
        $tran->seq_region_start - 5000, $tran->seq_region_end + 5000);
}

sub within_feature {
    my $vfoa    = shift;
    my $vf      = $vfoa->variation_feature;
    my $feat    = $vfoa->feature;
    
    return overlap(
        $vf->seq_region_start, 
        $vf->seq_region_end,
        $feat->seq_region_start, 
        $feat->seq_region_end
    );
}

sub within_transcript {
    my $tva = shift;
    return within_feature($tva);
}

sub _before_start {
    my ($vf, $feat, $dist) = @_;
    return ( ($vf->seq_region_end >= ($feat->seq_region_start - $dist)) and 
        ($vf->seq_region_end < $feat->seq_region_start) );
}

sub _after_end {
    my ($vf, $feat, $dist) = @_;
    return ( ($vf->seq_region_start <= ($feat->seq_region_end + $dist)) 
            and ($vf->seq_region_start > $feat->seq_region_end) );
}

sub _upstream {
    my ($vf, $feat, $dist) = @_;
    return $feat->strand == 1 ? 
        _before_start($vf, $feat, $dist) : 
        _after_end($vf, $feat, $dist);
}

sub _downstream {
    my ($vf, $feat, $dist) = @_;
    return $feat->strand == 1 ? 
        _after_end($vf, $feat, $dist) : 
        _before_start($vf, $feat, $dist);
}

sub upstream_5KB {
    my $vfo     = shift;
    my $vf      = $vfo->variation_feature;
    my $feat    = $vfo->feature;

    return (_upstream($vf, $feat, 5000) and not upstream_2KB($vfo));
}

sub downstream_5KB {
    my $vfo     = shift;
    my $vf      = $vfo->variation_feature;
    my $feat    = $vfo->feature;

    return (_downstream($vf, $feat, 5000) and not downstream_500B($vfo));
}

sub upstream_2KB {
    my $vfo     = shift;
    my $vf      = $vfo->variation_feature;
    my $feat    = $vfo->feature; 

    return _upstream($vf, $feat, 2000);
}

sub downstream_2KB {
    my $vfo     = shift;
    my $vf      = $vfo->variation_feature;
    my $feat    = $vfo->feature; 

    return _downstream($vf, $feat, 2000);
}

sub downstream_500B {
    my $vfo     = shift;
    my $vf      = $vfo->variation_feature;
    my $feat    = $vfo->feature;

    return _downstream($vf, $feat, 500);
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
                    if (overlap($vf->seq_region_start, $vf->seq_region_end, 
                            $coord->start, $coord->end)) {
                        return 1;
                    }
                }
            }
        }
    }
    
    return 0;
}

sub within_miRNA_target {
    my $tva = shift;
    # XXX: implement me!
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
    
    return overlap($vf->seq_region_start, $vf->seq_region_end, 
        $tran->coding_region_start, $tran->coding_region_end);
}

sub _before_coding {
    my ($vf, $tran) = @_;
    return 0 unless defined $tran->translation;
    return overlap($vf->seq_region_start, $vf->seq_region_end,
        $tran->seq_region_start, $tran->coding_region_start-1);    
}

sub _after_coding {
    my ($vf, $tran) = @_;
    return 0 unless defined $tran->translation;
    return overlap($vf->seq_region_start, $vf->seq_region_end, 
        $tran->coding_region_end+1, $tran->seq_region_end);    
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
    
    return 0 unless $vf->var_class =~ /^(in|del)/;

    return @{ $tva->transcript_variation->cds_coords } > 1;
}

sub synonymous_coding {
    my $tva = shift;
    my $tv  = $tva->transcript_variation;

    my $alt_aa = $tva->peptide;

    return 0 unless defined $alt_aa;
    
    my $ref_aa = $tv->reference_allele->peptide;

    return ( $alt_aa eq $ref_aa );
}

sub non_synonymous_coding {
    my $tva = shift;
    my $tv  = $tva->transcript_variation;
    
    return 0 if frameshift($tva);
    return 0 if stop_gained($tva);

    my $alt_aa = $tva->peptide;
    
    return 0 unless defined $alt_aa;

    my $ref_aa = $tv->reference_allele->peptide;
    
    return ( $alt_aa ne $ref_aa );
}

sub stop_gained {
    my $tva = shift;
    my $tv  = $tva->transcript_variation;

    my $alt_aa = $tva->peptide;

    return 0 unless defined $alt_aa;
    
    my $ref_aa = $tv->reference_allele->peptide;

    return ( ($alt_aa =~ /\*/) and ($ref_aa !~ /\*/) );
}

sub stop_lost {
    my $tva = shift;
    my $tv  = $tva->transcript_variation;

    my $alt_aa = $tva->peptide;

    return 0 unless defined $alt_aa;
    
    my $ref_aa = $tv->reference_allele->peptide;

    return ( ($alt_aa !~ /\*/) and ($ref_aa =~ /\*/) );
}

sub frameshift {
    my $tva = shift;

    return 0 if partial_codon($tva);
    return 0 if coding_unknown($tva);

    my $tv = $tva->transcript_variation;

    return 0 unless defined $tv->cds_start;
    
    my $var_len = $tv->cds_end - $tv->cds_start + 1;

    my $seq = $tva->feature_seq;
    
    $seq = '' if $seq eq '-';

    return abs( length($seq) - $var_len ) % 3;
}

sub partial_codon {
    my $tva = shift;
    
    my $tv = $tva->transcript_variation;
    
    return 0 unless defined $tv->translation_start;

    my $cds_length = length $tv->translateable_seq;

    my $codon_cds_start = ($tv->translation_start * 3) - 2;

    my $last_codon_length = $cds_length - ($codon_cds_start - 1);
    
    return ( $last_codon_length < 3 and $last_codon_length > 0 );
}

sub within_coding_frameshift_intron {
    my $tva = shift;
    
    return (within_coding_region($tva) and 
        $tva->transcript_variation->intron_effects->{within_frameshift_intron});
}

sub affects_peptide {
    my $tva     = shift;
    return (within_coding_region($tva) and (not within_intron($tva)));
}

sub coding_unknown {
    my $tva = shift;
    return (affects_peptide($tva) and ($tva->allele_string !~ /\//));
}

sub within_regulatory_feature {
    my $rfva = shift;
    return within_feature($rfva);
}

sub within_motif_feature {
    my $mfva = shift;
    return within_feature($mfva);
}

sub increased_binding_affinity {
    my $mfva = shift;
    return (within_motif_feature($mfva) and ($mfva->binding_affinity_change > 0));
}

sub decreased_binding_affinity {
    my $mfva = shift;
    return (within_motif_feature($mfva) and ($mfva->binding_affinity_change < 0));   
}

1;

