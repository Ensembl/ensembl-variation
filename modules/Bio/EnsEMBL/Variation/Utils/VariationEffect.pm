=head1 LICENSE

 Copyright (c) 1999-2012 The European Bioinformatics Institute and
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

=head1 NAME

Bio::EnsEMBL::Variation::Utils::VariationEffect

=head1 DESCRIPTION

This module defines a set of predicate subroutines that check the effect of a
Bio::EnsEMBL::Variation::VariationFeature on some other Bio::EnsEMBL::Feature. 
All of these predicates take a VariationFeatureOverlapAllele as their first and
only argument and return a true or false value depending on whether the effect
being checked for holds or not. The link between these predicates and the 
specific effect is configured in the Bio::EnsEMBL::Variation::Utils::Config 
module and a list of OverlapConsequence objects that represent a link between,
for example, a Sequence Ontology consequence term, and the predicate that
checks for it is provided in the Bio::EnsEMBL::Variation::Utils::Constants
module. If you want to add a new consequence you should write a predicate in 
this module and then add an entry in the configuration file.

=cut

package Bio::EnsEMBL::Variation::Utils::VariationEffect;

use strict;
use warnings;

use base qw(Exporter);

our @EXPORT_OK = qw(overlap within_cds MAX_DISTANCE_FROM_TRANSCRIPT);

use constant MAX_DISTANCE_FROM_TRANSCRIPT => 5000;

#package Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele;

sub overlap {
    my ( $f1_start, $f1_end, $f2_start, $f2_end ) = @_;
   
    return ( ($f1_end >= $f2_start) and ($f1_start <= $f2_end) );
}

sub within_feature {
    my $bvfoa   = shift;
    my $bvf     = $bvfoa->base_variation_feature;
    my $feat    = $bvfoa->feature;
    
    return overlap(
        $bvf->start, 
        $bvf->end,
        $feat->start, 
        $feat->end
    );
}

sub partial_overlap_feature {
    my $bvfoa   = shift;
    my $bvf     = $bvfoa->base_variation_feature;
    my $feat    = $bvfoa->feature;
    
    return (
        within_feature($bvfoa) and 
        (not complete_overlap_feature($bvfoa)) and
        (($bvf->end > $feat->end) or ($bvf->start < $feat->start))
    );
}

sub complete_within_feature {
    my $bvfoa   = shift;
    my $bvf     = $bvfoa->base_variation_feature;
    my $feat    = $bvfoa->feature;
    
    return (
        ($bvf->start >= $feat->start) and 
        ($bvf->end <= $feat->end)
    );
}

sub complete_overlap_feature {
    my $bvfoa   = shift;
    my $bvf     = $bvfoa->base_variation_feature;
    my $feat    = $bvfoa->feature;
    
    return ( 
        ($bvf->start <= $feat->start) and 
        ($bvf->end >= $feat->end) 
    );
}

sub deletion {
    my $bvfoa = shift;
   
    return (
        ($bvfoa->base_variation_feature->class_SO_term eq 'deletion') or
        ($bvfoa->base_variation_feature->class_SO_term =~ /deletion/i) or
        ($bvfoa->base_variation_feature->class_SO_term =~ /loss/i)
    );
}

sub insertion {
    my $bvfoa = shift;
    
    return (
        ($bvfoa->base_variation_feature->class_SO_term eq 'insertion') or
        ($bvfoa->base_variation_feature->class_SO_term =~ /insertion/i)
    );
}

sub duplication {
    my $bvfoa = shift;
    
    return (
        (
            ($bvfoa->base_variation_feature->class_SO_term eq 'duplication') or
            ($bvfoa->base_variation_feature->class_SO_term =~ /duplication/i)
        ) and
        (not tandem_duplication($bvfoa))
    );
}

sub tandem_duplication {
    my $bvfoa = shift;
    
    return (
        ($bvfoa->base_variation_feature->class_SO_term eq 'tandem_duplication') or
        ($bvfoa->base_variation_feature->class_SO_term =~ /tandem_duplication/i) 
    );
}

sub _before_start {
    my ($vf, $feat, $dist) = @_;
    
    return ( ($vf->end >= ($feat->start - $dist)) and 
        ($vf->end < $feat->start) );
}

sub _after_end {
    my ($vf, $feat, $dist) = @_;
    return ( ($vf->start <= ($feat->end + $dist)) 
            and ($vf->start > $feat->end) );
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

#package Bio::EnsEMBL::Variation::TranscriptVariationAllele;

sub upstream_5KB {
    my $vfoa    = shift;
    my $vf      = $vfoa->variation_feature;
    my $feat    = $vfoa->feature;

    return (_upstream($vf, $feat, 5000) and not upstream_2KB($vfoa));
}

sub downstream_5KB {
    my $vfoa    = shift;
    my $vf      = $vfoa->variation_feature;
    my $feat    = $vfoa->feature;

    return (_downstream($vf, $feat, 5000) and not downstream_500B($vfoa));
}

sub upstream_2KB {
    my $vfoa    = shift;
    my $vf      = $vfoa->variation_feature;
    my $feat    = $vfoa->feature; 

    return _upstream($vf, $feat, 2000);
}

sub downstream_2KB {
    my $vfoa    = shift;
    my $vf      = $vfoa->variation_feature;
    my $feat    = $vfoa->feature; 

    return _downstream($vf, $feat, 2000);
}

sub downstream_500B {
    my $vfoa    = shift;
    my $vf      = $vfoa->variation_feature;
    my $feat    = $vfoa->feature;

    return _downstream($vf, $feat, 500);

}

sub affects_transcript {
    my ($vf, $tran) = @_;
    
    return 0 unless $tran->isa('Bio::EnsEMBL::Transcript');
    
    return overlap(
        $vf->start, 
        $vf->end,
        $tran->start - 5000, 
        $tran->end + 5000
    );
}

sub within_transcript {
    my $tva = shift;
    return within_feature($tva);
}

sub within_nmd_transcript {
    my $tva     = shift;
    my $tran    = $tva->transcript; 

    return ( within_transcript($tva) and ($tran->biotype eq 'nonsense_mediated_decay') );
}

sub within_non_coding_gene {
    my $tva     = shift;
    my $tran    = $tva->transcript;
    
    return ( within_transcript($tva) and (not $tran->translation) and (not within_mature_miRNA($tva)));
}

sub within_miRNA {
    my $tva     = shift;
    my $tran    = $tva->transcript;
    
    # don't call this for now
    
    return 0;
    
    return ( within_transcript($tva) and ($tran->biotype eq 'miRNA') );
}

sub within_mature_miRNA {
    my $tva     = shift;
    
    my $tv      = $tva->transcript_variation;
    my $vf      = $tva->variation_feature;
    my $tran    = $tva->transcript;
        
    return 0 unless ( within_transcript($tva) and ($tran->biotype eq 'miRNA') );
        
    my ($attribute) = @{ $tran->get_all_Attributes('miRNA') };
    
    if (defined $attribute && $attribute->value =~ /(\d+)-(\d+)/) { 
        for my $coord ($tv->_mapper->cdna2genomic($1, $2, $tran->strand)) {
            if ($coord->isa('Bio::EnsEMBL::Mapper::Coordinate')) {
                if (overlap(
                        $vf->start, 
                        $vf->end, 
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
        $tva->transcript_variation->_intron_effects->{start_splice_site} :
        $tva->transcript_variation->_intron_effects->{end_splice_site};
}

sub acceptor_splice_site {
    my $tva     = shift;
    my $tran    = $tva->transcript;
    
    return $tran->strand == 1 ? 
        $tva->transcript_variation->_intron_effects->{end_splice_site} :
        $tva->transcript_variation->_intron_effects->{start_splice_site};
}

sub essential_splice_site {
    my $tva = shift;
    
    return ( acceptor_splice_site($tva) or donor_splice_site($tva) );
}

sub splice_region {
    my $tva    = shift;

    return $tva->transcript_variation->_intron_effects->{splice_region};
}

sub within_intron {
    my $tva    = shift;

    return $tva->transcript_variation->_intron_effects->{intronic};
}

sub within_cds {
    my $tva     = shift;
    my $vf      = $tva->variation_feature;
    my $tran    = $tva->transcript;
    my $tv      = $tva->transcript_variation;
    
    my $cds_coords = $tva->transcript_variation->cds_coords;
    
    if (@$cds_coords > 0) {
        for my $coord (@$cds_coords) {
            if ($coord->isa('Bio::EnsEMBL::Mapper::Coordinate')) {
                if ($coord->end > 0 && $coord->start <= length($tv->_translateable_seq)) { 
                    return 1;
                }
            }
        }
    }

    # we also need to check if the vf is in a frameshift intron within the CDS

    if (defined $tran->translation &&
        $tva->transcript_variation->_intron_effects->{within_frameshift_intron}) {
 
        return overlap(
            $vf->start, 
            $vf->end, 
            $tran->coding_region_start,
            $tran->coding_region_end,
        );
    }
        
    return 0;
}

sub within_cdna {
    my $tva     = shift;
    my $vf      = $tva->variation_feature;
    my $tran    = $tva->transcript;
    my $tv      = $tva->transcript_variation;
    
    my $cdna_coords = $tv->cdna_coords;
    
    if (@$cdna_coords > 0) {
        for my $coord (@$cdna_coords) {
            if ($coord->isa('Bio::EnsEMBL::Mapper::Coordinate')) {
                if ($coord->end > 0 && $coord->start <= $tran->length) {
                    return 1;
                }
            }
        }
    }
    
    # we also need to check if the vf is in a frameshift intron within the cDNA

    if ($tva->transcript_variation->_intron_effects->{within_frameshift_intron}) {
        return within_transcript($tva); 
    }
    
    return 0;
}

sub _before_coding {
    my ($vf, $tran) = @_;
    return 0 unless defined $tran->translation;
    
    my $vf_s  = $vf->start;
    my $vf_e  = $vf->end;
    my $t_s   = $tran->start;
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
    
    my $vf_s  = $vf->start;
    my $vf_e  = $vf->end;
    my $t_e   = $tran->end;
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
   
    # pass the no_db flag to var_class to ensure we don't rely on the database for it 
    # as it may not have been set at this stage in the pipeline
    my $class = $vf->var_class(1);

    return 0 unless $class =~ /insertion|deletion|indel/;

    return @{ $tva->transcript_variation->cds_coords } > 1;
}

sub _get_peptide_alleles {
    my $tva = shift;
    my $tv  = $tva->transcript_variation;
    
    return () if frameshift($tva);

    my $alt_pep = $tva->peptide;
    
    return () unless defined $alt_pep;
    
    my $ref_pep = $tv->get_reference_TranscriptVariationAllele->peptide;
    
    $ref_pep = '' if $ref_pep eq '-';
    $alt_pep = '' if $alt_pep eq '-';
    
    return ($ref_pep, $alt_pep);
}

sub _get_codon_alleles {
    my $tva = shift;
    my $tv  = $tva->transcript_variation;
    
    return () if frameshift($tva);

    my $alt_codon = $tva->codon;
    
    return () unless defined $alt_codon;
    
    my $ref_codon = $tv->get_reference_TranscriptVariationAllele->codon;
    
    $ref_codon = '' if $ref_codon eq '-';
    $alt_codon = '' if $alt_codon eq '-';
    
    return ($ref_codon, $alt_codon);
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

    return ( ($alt_pep eq $ref_pep) and (not stop_retained($tva) and ($alt_pep !~ /X/) and ($ref_pep !~ /X/)) );
}

sub non_synonymous_codon {
    my $tva = shift;
    
    my ($ref_pep, $alt_pep) = _get_peptide_alleles($tva);
    
    return 0 unless defined $ref_pep;
    
    return 0 if affects_start_codon($tva);
    return 0 if stop_lost($tva);
    return 0 if stop_gained($tva);
    return 0 if partial_codon($tva);

    return 0 if inframe_codon_loss($tva);
    return 0 if inframe_codon_gain($tva);
    
    return ( $ref_pep ne $alt_pep );
}

sub inframe_codon_gain {
    my $tva = shift;
    
    my ($ref_codon, $alt_codon) = _get_codon_alleles($tva);
    
    return 0 unless defined $ref_codon;
    
    return ( 
        (length($alt_codon) > length ($ref_codon)) &&
        ( ($alt_codon =~ /^\Q$ref_codon\E/) || ($alt_codon =~ /\Q$ref_codon\E$/) )
    );
}

sub inframe_codon_loss {
    my $tva = shift;
    
    my ($ref_codon, $alt_codon) = _get_codon_alleles($tva);

    return 0 unless defined $ref_codon;
  
    return ( 
        (length($alt_codon) < length ($ref_codon)) &&
        ( ($ref_codon =~ /^\Q$alt_codon\E/) || ($ref_codon =~ /\Q$alt_codon\E$/) )
    );
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

    my $allele_len = $tva->seq_length;

    # if the allele length is undefined then we can't call a frameshift

    return 0 unless defined $allele_len;

    return abs( $allele_len - $var_len ) % 3;
}

sub partial_codon {
    my $tva = shift;
    
    my $tv = $tva->transcript_variation;
    
    return 0 unless defined $tv->translation_start;

    my $cds_length = length $tv->_translateable_seq;

    my $codon_cds_start = ($tv->translation_start * 3) - 2;

    my $last_codon_length = $cds_length - ($codon_cds_start - 1);
    
    return ( $last_codon_length < 3 and $last_codon_length > 0 );
}

sub coding_unknown {
    my $tva = shift;

    return (within_cds($tva) and ((not $tva->peptide) or ($tva->peptide =~ /X/)) and (not frameshift($tva)));
}

#package Bio::EnsEMBL::Variation::RegulatoryFeatureVariationAllele;

sub within_regulatory_feature {
    my $rfva = shift;
    return within_feature($rfva);
}

#package Bio::EnsEMBL::Variation::ExternalFeatureVariationAllele;

sub within_external_feature {
    my $efva = shift;
    return (within_feature($efva) and (not within_miRNA_target_site($efva)));
}

#sub within_miRNA_target_site {
#    my $efva = shift;
#
#    my $fset = $efva->variation_feature_overlap->feature->feature_set;
#
#    return ($fset && $fset->name eq 'miRanda miRNA targets');
#}

#package Bio::EnsEMBL::Variation::MotifFeatureVariationAllele;

#sub within_motif_feature {
#    my $mfva = shift;
#    return (
#        within_feature($mfva) and
#        !increased_binding_affinity($mfva) and
#        !decreased_binding_affinity($mfva) 
#    );
#}

sub within_motif_feature {
    my $mfva = shift;
    return within_feature($mfva);
}

#sub increased_binding_affinity {
#    my $mfva = shift;
#    my $change = $mfva->binding_affinity_change;
#    return (within_feature($mfva) and (defined $change) and ($change > 0));
#}
#
#sub decreased_binding_affinity {
#    my $mfva = shift;
#    my $change = $mfva->binding_affinity_change;
#    return (within_feature($mfva) and (defined $change) and ($change < 0));
#}

sub contains_entire_feature {
    my $vfo = shift;

    my $vf = $vfo->variation_feature;
    my $feat = $vfo->feature;

    return ( ($vf->start <= $feat->start) && ($vf->end >= $feat->end) ); 
}

1;

