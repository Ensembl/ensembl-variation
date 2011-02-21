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

our @EXPORT_OK = qw(overlap within_cds);

sub overlap {
    my ( $f1_start, $f1_end, $f2_start, $f2_end ) = @_;
   
    return ( ($f1_end >= $f2_start) and ($f1_start <= $f2_end) );
}

sub affects_transcript {
    my ($vf, $tran) = @_;
    
    return 0 unless $tran->isa('Bio::EnsEMBL::Transcript');
    
    return overlap(
        $vf->seq_region_start, 
        $vf->seq_region_end,
        $tran->seq_region_start - 5000, 
        $tran->seq_region_end + 5000
    );
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
    my $tran    = $tva->transcript;
    
    return ( within_transcript($tva) and ($tran->biotype eq 'miRNA') );
}

sub within_mature_miRNA {
    my $tva     = shift;
    
    return 0 unless within_miRNA($tva);
    
    my $tv      = $tva->transcript_variation;
    my $vf      = $tva->variation_feature;
    my $tran    = $tva->transcript;
        
    my ($attribute) = @{ $tran->get_all_Attributes('miRNA') };
    
    if (defined $attribute && $attribute->value =~ /(\d+)-(\d+)/) { 
        for my $coord ($tv->mapper->cdna2genomic($1, $2, $tran->strand)) {
            if ($coord->isa('Bio::EnsEMBL::Mapper::Coordinate')) {
                if (overlap(
                        $vf->seq_region_start, 
                        $vf->seq_region_end, 
                        $coord->start, 
                        $coord->end) ) {
                    return 1;
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

sub within_cds {
    my $tva     = shift;
    my $vf      = $tva->variation_feature;
    my $tran    = $tva->transcript;
    
    my $cds_coords = $tva->transcript_variation->cds_coords;
    
    return ( 
        (@$cds_coords > 0) and 
        ( grep {$_->isa('Bio::EnsEMBL::Mapper::Coordinate')} @$cds_coords )
    );
}

sub within_cdna {
    my $tva     = shift;
    my $vf      = $tva->variation_feature;
    my $tran    = $tva->transcript;
    
    my $cdna_coords = $tva->transcript_variation->cdna_coords;
    
    return ( 
        (@$cdna_coords > 0) and 
        ( grep {$_->isa('Bio::EnsEMBL::Mapper::Coordinate')} @$cdna_coords )
    );
}

sub _before_coding {
    my ($vf, $tran) = @_;
    return 0 unless defined $tran->translation;
    
    my $vf_s  = $vf->seq_region_start;
    my $vf_e  = $vf->seq_region_end;
    my $t_s   = $tran->seq_region_start;
    my $cds_s = $tran->coding_region_start;
    
    # we need to special case insertions just before the CDS start
    if ($vf_s == $vf_e+1 && $vf_s == $cds_s) {
        return 1;
    }
    
    return overlap($vf_s, $vf_e, $t_s, $cds_s-1);    
}

sub _after_coding {
    my ($vf, $tran) = @_;
    return 0 unless defined $tran->translation;
    
    my $vf_s  = $vf->seq_region_start;
    my $vf_e  = $vf->seq_region_end;
    my $t_e   = $tran->seq_region_end;
    my $cds_e = $tran->coding_region_end;
    
    # we need to special case insertions just after the CDS end
    if ($vf_s == $vf_e+1 && $vf_e == $cds_e) {
        return 1;
    }
    
    return overlap($vf_s, $vf_e, $cds_e+1, $t_e);
}

sub within_5_prime_utr {
    my $tva     = shift;
    my $vf      = $tva->variation_feature;
    my $tran    = $tva->transcript;
    
    my $five_prime_of_coding = 
        $tran->strand == 1 ? 
            _before_coding($vf, $tran) : 
            _after_coding($vf, $tran);
    
    return ( $five_prime_of_coding and within_cdna($tva) );
}

sub within_3_prime_utr {
    my $tva     = shift;
    my $vf      = $tva->variation_feature;
    my $tran    = $tva->transcript; 
    
    my $three_prime_of_coding = 
        $tran->strand == 1 ? 
            _after_coding($vf, $tran) : 
            _before_coding($vf, $tran);
    
    return ( $three_prime_of_coding and within_cdna($tva) );
}

sub complex_indel {
    my $tva     = shift;
    my $vf      = $tva->variation_feature;
    
    return 0 unless $vf->var_class =~ /^(in|del)/;

    return @{ $tva->transcript_variation->cds_coords } > 1;
}

sub _get_peptide_alleles {
    my $tva = shift;
    my $tv  = $tva->transcript_variation;
    
    return () if frameshift($tva);

    my $alt_pep = $tva->peptide;
    
    return () unless defined $alt_pep;
    
    my $ref_pep = $tv->reference_allele->peptide;
    
    $ref_pep = '' if $ref_pep eq '-';
    $alt_pep = '' if $alt_pep eq '-';
    
    return ($ref_pep, $alt_pep);
}

sub stop_retained {
    my $tva = shift;
    
    my ($ref_pep, $alt_pep) = _get_peptide_alleles($tva);
    
    return 0 unless $ref_pep;

    return ( $alt_pep =~ /\*/ && $ref_pep =~ /\*/ );
}

sub affects_start_codon {
    my $tva = shift;
    my $tv  = $tva->transcript_variation;
    
    my ($ref_pep, $alt_pep) = _get_peptide_alleles($tva);

    return 0 unless $ref_pep;

    return ( ($tv->translation_start == 1) and (substr($ref_pep,0,1) ne substr($alt_pep,0,1)) );
}

sub synonymous_codon {
    my $tva = shift;

    my ($ref_pep, $alt_pep) = _get_peptide_alleles($tva);
    
    return 0 unless $ref_pep;

    return ( ($alt_pep eq $ref_pep) and (not stop_retained($tva)) );
}

sub non_synonymous_codon {
    my $tva = shift;
    
    my ($ref_pep, $alt_pep) = _get_peptide_alleles($tva);
    
    return 0 unless defined $ref_pep;
    
    return 0 if affects_start_codon($tva);
    return 0 if stop_lost($tva);
    return 0 if partial_codon($tva);
    
    return ( ($ref_pep ne $alt_pep) and (length($ref_pep) == length($alt_pep)) );
}

sub inframe_codon_gain {
    my $tva = shift;
    
    my ($ref_pep, $alt_pep) = _get_peptide_alleles($tva);
    
    return 0 unless defined $ref_pep;
    
    return (length($ref_pep) < length($alt_pep));
}

sub inframe_codon_loss {
    my $tva = shift;
    
    my ($ref_pep, $alt_pep) = _get_peptide_alleles($tva);
    
    return 0 unless defined $ref_pep;
    
    return (length($ref_pep) > length($alt_pep));
}

sub stop_gained {
    my $tva = shift;
    my $tv  = $tva->transcript_variation;

    my ($ref_pep, $alt_pep) = _get_peptide_alleles($tva);
    
    return 0 unless defined $ref_pep;

    return ( ($alt_pep =~ /\*/) and ($ref_pep !~ /\*/) );
}

sub stop_lost {
    my $tva = shift;

    my ($ref_pep, $alt_pep) = _get_peptide_alleles($tva);
    
    return 0 unless defined $ref_pep;

    return ( ($alt_pep !~ /\*/) and ($ref_pep =~ /\*/) );
}

sub frameshift {
    my $tva = shift;

    return 0 if partial_codon($tva);

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
    
    return (within_cds($tva) and 
        $tva->transcript_variation->intron_effects->{within_frameshift_intron});
}

sub coding_unknown {
    my $tva = shift;
      
    return (within_cds($tva) and (not $tva->peptide) and (not frameshift($tva)));
}


sub within_miRNA_target_site {
    my $tva = shift;
    # XXX: implement me!
    return 0;
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

