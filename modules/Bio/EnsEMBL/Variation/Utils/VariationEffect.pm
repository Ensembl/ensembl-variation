=head1 LICENSE

Copyright 2013 Ensembl

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut


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

our @EXPORT_OK = qw(overlap within_cds MAX_DISTANCE_FROM_TRANSCRIPT within_intron stop_lost affects_start_codon $UPSTREAM_DISTANCE $DOWNSTREAM_DISTANCE);

use constant MAX_DISTANCE_FROM_TRANSCRIPT => 5000;

our $UPSTREAM_DISTANCE = MAX_DISTANCE_FROM_TRANSCRIPT;
our $DOWNSTREAM_DISTANCE = MAX_DISTANCE_FROM_TRANSCRIPT;

#package Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele;

sub overlap {
    my ( $f1_start, $f1_end, $f2_start, $f2_end ) = @_;
   
    return ( ($f1_end >= $f2_start) and ($f1_start <= $f2_end) );
}

sub within_feature {
    my ($bvfoa, $feat) = @_;
    my $bvf            = $bvfoa->base_variation_feature;
    $feat              = $bvfoa->feature unless defined($feat);
    
    return overlap(
        $bvf->start, 
        $bvf->end,
        $feat->start, 
        $feat->end
    );
}

sub partial_overlap_feature {
    my ($bvfoa, $feat) = @_;
    my $bvf            = $bvfoa->base_variation_feature;
    $feat              = $bvfoa->feature unless defined($feat);
    
    return (
        within_feature(@_) and 
        (not complete_overlap_feature(@_)) and
        (($bvf->end > $feat->end) or ($bvf->start < $feat->start))
    );
}

sub complete_within_feature {
    my ($bvfoa, $feat) = @_;
    my $bvf            = $bvfoa->base_variation_feature;
    $feat              = $bvfoa->feature unless defined($feat);
    
    return (
        ($bvf->start >= $feat->start) and 
        ($bvf->end <= $feat->end)
    );
}

sub complete_overlap_feature {
    my ($bvfoa, $feat) = @_;
    my $bvf            = $bvfoa->base_variation_feature;
    $feat              = $bvfoa->feature unless defined($feat);
    
    return ( 
        ($bvf->start <= $feat->start) and 
        ($bvf->end >= $feat->end) 
    );
}

sub deletion {
    my $bvfoa = shift;
    
    my $bvf = $bvfoa->base_variation_feature;
    
    # sequence variant will have alleles
    if($bvf->isa('Bio::EnsEMBL::Variation::VariationFeature')) {
        my ($ref_allele, $alt_allele) = _get_alleles($bvfoa);
        return (
            (defined($ref_allele) && ($alt_allele eq '') and $ref_allele) or
            $bvf->allele_string =~ /deletion/i
        );
    }
    
    # structural variant depends on class
    if($bvf->isa('Bio::EnsEMBL::Variation::StructuralVariationFeature')) {
        return (
            ($bvfoa->base_variation_feature->class_SO_term(undef, 1) eq 'deletion') or
            ($bvfoa->base_variation_feature->class_SO_term(undef, 1) =~ /deletion/i) or
            ($bvfoa->base_variation_feature->class_SO_term(undef, 1) =~ /loss/i)
        );
    }
    
    else { return 0; }
}

sub insertion {
    my $bvfoa = shift;
    
    my $bvf = $bvfoa->base_variation_feature;
    
    # sequence variant will have alleles
    if($bvf->isa('Bio::EnsEMBL::Variation::VariationFeature')) {
        my ($ref_allele, $alt_allele) = _get_alleles($bvfoa);
        return (
            (defined($ref_allele) && ($ref_allele eq '') and $alt_allele) or
            $bvf->allele_string =~ /insertion/i
        );
    }
    
    # structural variant depends on class
    if($bvf->isa('Bio::EnsEMBL::Variation::StructuralVariationFeature')) {
        return (
            duplication($bvfoa) or
            tandem_duplication($bvfoa) or
            ($bvfoa->base_variation_feature->class_SO_term(undef, 1) eq 'insertion') or
            ($bvfoa->base_variation_feature->class_SO_term(undef, 1) =~ /insertion/i) or
            ($bvfoa->base_variation_feature->class_SO_term(undef, 1) =~ /gain/i)
        );
    }
    
    else { return 0; }
}

sub copy_number_gain {
    my $bvfoa = shift;
    
    return (duplication($bvfoa) or tandem_duplication($bvfoa) or $bvfoa->base_variation_feature->class_SO_term(undef, 1) =~ /gain/i);
}

sub copy_number_loss {
    my $bvfoa = shift;
    
    return $bvfoa->base_variation_feature->class_SO_term(undef, 1) =~ /loss/i;
}

sub duplication {
    my $bvfoa = shift;
    
    return (
        (
            ($bvfoa->base_variation_feature->class_SO_term(undef, 1) eq 'duplication') or
            ($bvfoa->base_variation_feature->class_SO_term(undef, 1) =~ /duplication/i)
        ) and
        (not tandem_duplication($bvfoa))
    );
}

sub tandem_duplication {
    my $bvfoa = shift;
    
    my $bvf = $bvfoa->base_variation_feature;
    
    # for sequence variants, check sequence vs ref
    if($bvf->isa('Bio::EnsEMBL::Variation::VariationFeature')) {
        my ($ref_allele, $alt_allele) = _get_alleles($bvfoa);
        
        return 0 unless $ref_allele and $alt_allele;
        return 0 unless
            length($alt_allele) > length($ref_allele) and
            length($alt_allele) % length($ref_allele) == 0;
        
        my $copies = length($alt_allele) / length($ref_allele);
        
        return $alt_allele eq $ref_allele x $copies;
    }
    
    # structural variant depends on class
    if($bvf->isa('Bio::EnsEMBL::Variation::StructuralVariationFeature')) {
        return (
            ($bvfoa->base_variation_feature->class_SO_term(undef, 1) eq 'tandem_duplication') or
            ($bvfoa->base_variation_feature->class_SO_term(undef, 1) =~ /tandem_duplication/i) 
        );
    }
}

sub feature_ablation {
    my $bvfoa = shift;
    
    return (deletion($bvfoa) and complete_overlap_feature($bvfoa));
}

sub feature_amplification {
    my $bvfoa = shift;
    
    return (copy_number_gain($bvfoa) && complete_overlap_feature($bvfoa));
}

sub feature_elongation {
    my $bvfoa = shift;
    
    return (
        complete_within_feature($bvfoa) and
        (copy_number_gain($bvfoa) or insertion($bvfoa)) and
        not(
            $bvfoa->isa('Bio::EnsEMBL::Variation::BaseTranscriptVariationAllele') and
            (inframe_insertion($bvfoa) or stop_lost($bvfoa))
        )
    );
}

sub feature_truncation {
    my $bvfoa = shift;
    
    return (
        (partial_overlap_feature($bvfoa) or complete_within_feature($bvfoa)) and
        (copy_number_loss($bvfoa) or deletion($bvfoa)) and
        not(
            $bvfoa->isa('Bio::EnsEMBL::Variation::BaseTranscriptVariationAllele') and
            (inframe_deletion($bvfoa) or stop_gained($bvfoa))
        )
    );
}

#sub transcript_fusion {
#    #my $bvfoa = shift;
#    #my $bvf   = $bvfoa->base_variation_feature;
#    
#    return 0;
#    
#    #my $transcripts = $bvf->_get_overlapping_Transcripts();
#}

sub _before_start {
    my ($bvf, $feat, $dist) = @_;
    
    return ( ($bvf->end >= ($feat->start - $dist)) and 
        ($bvf->end < $feat->start) );
}

sub _after_end {
    my ($bvf, $feat, $dist) = @_;
    return ( ($bvf->start <= ($feat->end + $dist)) 
            and ($bvf->start > $feat->end) );
}

sub _upstream {
    my ($bvf, $feat, $dist) = @_;
    return $feat->strand == 1 ? 
        _before_start($bvf, $feat, $dist) : 
        _after_end($bvf, $feat, $dist);
}

sub _downstream {
    my ($bvf, $feat, $dist) = @_;
    return $feat->strand == 1 ? 
        _after_end($bvf, $feat, $dist) : 
        _before_start($bvf, $feat, $dist);
}

#package Bio::EnsEMBL::Variation::TranscriptVariationAllele;

sub upstream {
    my $vfoa    = shift;
    my $bvf     = $vfoa->base_variation_feature;
    my $feat    = $vfoa->feature;

    return _upstream($bvf, $feat, $UPSTREAM_DISTANCE);
}

sub downstream {
    my $vfoa    = shift;
    my $bvf     = $vfoa->base_variation_feature;
    my $feat    = $vfoa->feature;

    return _downstream($bvf, $feat, $DOWNSTREAM_DISTANCE);
}

sub affects_transcript {
    my ($bvf, $tran) = @_;
    
    return 0 unless $tran->isa('Bio::EnsEMBL::Transcript');
    
    return overlap(
        $bvf->start, 
        $bvf->end,
        $tran->start - 5000, 
        $tran->end + 5000
    );
}

sub within_transcript {
    my $bvfoa = shift;
    return within_feature($bvfoa);
}

sub within_nmd_transcript {
    my $bvfoa   = shift;
    my $tran    = $bvfoa->transcript; 

    return ( within_transcript($bvfoa) and ($tran->biotype eq 'nonsense_mediated_decay') );
}

sub within_non_coding_gene {
    my $bvfoa   = shift;
    my $tran    = $bvfoa->transcript;
    
    return ( within_transcript($bvfoa) and (not $tran->translation) and (not within_mature_miRNA($bvfoa)));
}

sub non_coding_exon_variant {
    my $bvfoa = shift;
    
    return 0 unless within_non_coding_gene($bvfoa);
    
    my $bvf   = $bvfoa->base_variation_feature;
    my $exons = $bvfoa->base_variation_feature_overlap->_exons;
    
    if(scalar grep {overlap($bvf->start, $bvf->end, $_->start, $_->end)} @$exons) {
        return 1;
    }
    else {
        return 0;
    }
}

sub within_miRNA {
    my $bvfoa   = shift;
    my $tran    = $bvfoa->transcript;
    
    # don't call this for now
    
    return 0;
    
    return ( within_transcript($bvfoa) and ($tran->biotype eq 'miRNA') );
}

sub within_mature_miRNA {
    my $bvfoa     = shift;
    
    my $bvfo      = $bvfoa->base_variation_feature_overlap;
    my $bvf       = $bvfoa->base_variation_feature;
    my $tran      = $bvfoa->transcript;
        
    return 0 unless ( within_transcript($bvfoa) and ($tran->biotype eq 'miRNA') );
        
    my ($attribute) = @{ $tran->get_all_Attributes('miRNA') };
    
    if (defined $attribute && $attribute->value =~ /(\d+)-(\d+)/) { 
        for my $coord ($bvfo->_mapper->cdna2genomic($1, $2, $tran->strand)) {
            if ($coord->isa('Bio::EnsEMBL::Mapper::Coordinate')) {
                if (overlap(
                        $bvf->start, 
                        $bvf->end, 
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
    my $bvfoa     = shift;
    my $tran      = $bvfoa->transcript;
    
    return $tran->strand == 1 ? 
        $bvfoa->base_variation_feature_overlap->_intron_effects->{start_splice_site} :
        $bvfoa->base_variation_feature_overlap->_intron_effects->{end_splice_site};
}

sub acceptor_splice_site {
    my $bvfoa     = shift;
    my $tran      = $bvfoa->transcript;
    
    return $tran->strand == 1 ? 
        $bvfoa->base_variation_feature_overlap->_intron_effects->{end_splice_site} :
        $bvfoa->base_variation_feature_overlap->_intron_effects->{start_splice_site};
}

sub essential_splice_site {
    my $bvfoa = shift;
    
    return ( acceptor_splice_site($bvfoa) or donor_splice_site($bvfoa) );
}

sub splice_region {
    my $bvfoa    = shift;
    
    return 0 if donor_splice_site($bvfoa);
    return 0 if acceptor_splice_site($bvfoa);
    return 0 if essential_splice_site($bvfoa);

    return $bvfoa->base_variation_feature_overlap->_intron_effects->{splice_region};
}

sub within_intron {
    my $bvfoa    = shift;

    return $bvfoa->base_variation_feature_overlap->_intron_effects->{intronic};
}

sub within_cds {
    my $bvfoa     = shift;
    my $bvf       = $bvfoa->base_variation_feature;
    my $tran      = $bvfoa->transcript;
    my $bvfo      = $bvfoa->base_variation_feature_overlap;
    
    my $cds_coords = $bvfoa->base_variation_feature_overlap->cds_coords;
    
    if (@$cds_coords > 0) {
        for my $coord (@$cds_coords) {
            if ($coord->isa('Bio::EnsEMBL::Mapper::Coordinate')) {
                if ($coord->end > 0 && $coord->start <= length($bvfo->_translateable_seq)) { 
                    return 1;
                }
            }
        }
    }

    # we also need to check if the vf is in a frameshift intron within the CDS

    if (defined $tran->translation &&
        $bvfoa->base_variation_feature_overlap->_intron_effects->{within_frameshift_intron}) {
 
        return overlap(
            $bvf->start, 
            $bvf->end, 
            $tran->coding_region_start,
            $tran->coding_region_end,
        );
    }
        
    return 0;
}

sub within_cdna {
    my $bvfoa     = shift;
    my $bvf       = $bvfoa->base_variation_feature;
    my $tran      = $bvfoa->transcript;
    my $bvfo      = $bvfoa->base_variation_feature_overlap;
    
    my $cdna_coords = $bvfo->cdna_coords;
    
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

    if ($bvfoa->base_variation_feature_overlap->_intron_effects->{within_frameshift_intron}) {
        return within_transcript($bvfoa); 
    }
    
    return 0;
}

sub _before_coding {
    my ($bvf, $tran) = @_;
    return 0 unless defined $tran->translation;
    
    my $bvf_s  = $bvf->start;
    my $bvf_e  = $bvf->end;
    my $t_s    = $tran->start;
    my $cds_s  = $tran->coding_region_start;
    
    # we need to special case insertions just before the CDS start
    if ($bvf_s == $bvf_e+1 && $bvf_s == $cds_s) {
        return 1;
    }
   
    return overlap($bvf_s, $bvf_e, $t_s, $cds_s-1);    
}

sub _after_coding {
    my ($bvf, $tran) = @_;
    return 0 unless defined $tran->translation;
    
    my $bvf_s  = $bvf->start;
    my $bvf_e  = $bvf->end;
    my $t_e    = $tran->end;
    my $cds_e  = $tran->coding_region_end;
    
    # we need to special case insertions just after the CDS end
    if ($bvf_s == $bvf_e+1 && $bvf_e == $cds_e) {
        return 1;
    }
    
    return overlap($bvf_s, $bvf_e, $cds_e+1, $t_e);
}

sub within_5_prime_utr {
    my $bvfoa     = shift;
    my $bvf       = $bvfoa->base_variation_feature;
    my $tran      = $bvfoa->transcript;
    
    my $five_prime_of_coding = 
        $tran->strand == 1 ? 
            _before_coding($bvf, $tran) : 
            _after_coding($bvf, $tran);
    
    return ( $five_prime_of_coding and within_cdna($bvfoa) );
}

sub within_3_prime_utr {
    my $bvfoa     = shift;
    my $bvf       = $bvfoa->base_variation_feature;
    my $tran      = $bvfoa->transcript; 
    
    my $three_prime_of_coding = 
        $tran->strand == 1 ? 
            _after_coding($bvf, $tran) : 
            _before_coding($bvf, $tran);
    
    return ( $three_prime_of_coding and within_cdna($bvfoa) );
}

sub complex_indel {
    my $bvfoa     = shift;
    my $bvf       = $bvfoa->base_variation_feature;
   
    # pass the no_db flag to var_class to ensure we don't rely on the database for it 
    # as it may not have been set at this stage in the pipeline
    my $class = $bvf->var_class(1);

    return 0 unless $class =~ /insertion|deletion|indel/;

    return @{ $bvfoa->base_variation_feature_overlap->cds_coords } > 1;
}

sub _get_peptide_alleles {
    my $bvfoa = shift;
    my $bvfo  = $bvfoa->base_variation_feature_overlap;
    
    return () if frameshift($bvfoa);

    my $alt_pep = $bvfoa->peptide;
    
    return () unless defined $alt_pep;
    
    my $ref_pep = $bvfo->get_reference_TranscriptVariationAllele->peptide;
    
    return () unless defined $ref_pep;
    
    $ref_pep = '' if $ref_pep eq '-';
    $alt_pep = '' if $alt_pep eq '-';
    
    return ($ref_pep, $alt_pep);
}

sub _get_codon_alleles {
    my $bvfoa = shift;
    my $bvfo  = $bvfoa->base_variation_feature_overlap;
    
    return () if frameshift($bvfoa);

    my $alt_codon = $bvfoa->codon;
    
    return () unless defined $alt_codon;
    
    my $ref_codon = $bvfo->get_reference_TranscriptVariationAllele->codon;
    
    return () unless defined $ref_codon;
    
    $ref_codon = '' if $ref_codon eq '-';
    $alt_codon = '' if $alt_codon eq '-';
    
    return ($ref_codon, $alt_codon);
}

sub _get_alleles {
    my $bvfoa = shift;
    my $bvfo  = $bvfoa->base_variation_feature_overlap;
    
    my $ref_tva = $bvfo->get_reference_VariationFeatureOverlapAllele;
    
    return () unless defined ($ref_tva);
    
    my $ref_allele = $ref_tva->variation_feature_seq;
    my $alt_allele = $bvfoa->variation_feature_seq;
    
    return () unless defined($ref_allele) and defined($alt_allele);
    
    $ref_allele = '' if $ref_allele eq '-';
    $alt_allele = '' if $alt_allele eq '-';
    
    return ($ref_allele, $alt_allele);
}

sub stop_retained {
    my $bvfoa = shift;
    
    my ($ref_pep, $alt_pep) = _get_peptide_alleles($bvfoa);
    
    return 0 unless $ref_pep;

    return ( $alt_pep =~ /\*/ && $ref_pep =~ /\*/ );
}

sub affects_start_codon {
    my $bvfoa = shift;
    my $bvfo  = $bvfoa->base_variation_feature_overlap;
    
    # sequence variant
    if($bvfo->isa('Bio::EnsEMBL::Variation::TranscriptVariation')) {
        my ($ref_pep, $alt_pep) = _get_peptide_alleles($bvfoa);
    
        return 0 unless $ref_pep;
    
        return ( ($bvfo->translation_start == 1) and (substr($ref_pep,0,1) ne substr($alt_pep,0,1)) );
    }
    
    # structural variant
    elsif($bvfo->isa('Bio::EnsEMBL::Variation::TranscriptStructuralVariation')) {
        my $tr  = $bvfo->transcript;
        my $bvf = $bvfo->base_variation_feature;
        
        my ($tr_crs, $tr_cre) = ($tr->coding_region_start, $tr->coding_region_end);
        return 0 unless defined($tr_crs) && defined($tr_cre);
        
        if($tr->strand == 1) {
            return overlap($tr_crs, $tr_crs + 2, $bvf->start, $bvf->end);
        }
        else {
            return overlap($tr_cre-2, $tr_cre, $bvf->start, $bvf->end);
        }
    }
    
    else {
        return 0;
    }
}

sub synonymous_variant {
    my $bvfoa = shift;

    my ($ref_pep, $alt_pep) = _get_peptide_alleles($bvfoa);
    
    return 0 unless $ref_pep;

    return ( ($alt_pep eq $ref_pep) and (not stop_retained($bvfoa) and ($alt_pep !~ /X/) and ($ref_pep !~ /X/)) );
}

sub missense_variant {
    my $bvfoa = shift;
    
    my ($ref_pep, $alt_pep) = _get_peptide_alleles($bvfoa);
    
    return 0 unless defined $ref_pep;
    
    return 0 if affects_start_codon($bvfoa);
    return 0 if stop_lost($bvfoa);
    return 0 if stop_gained($bvfoa);
    return 0 if partial_codon($bvfoa);

    return 0 if inframe_deletion($bvfoa);
    return 0 if inframe_insertion($bvfoa);
    
    return ( $ref_pep ne $alt_pep ) && ( length($ref_pep) == length($alt_pep) );
}

sub inframe_insertion {
    my $bvfoa = shift;
    
    # sequence variant
    if($bvfoa->base_variation_feature->isa('Bio::EnsEMBL::Variation::VariationFeature')) {
        my ($ref_codon, $alt_codon) = _get_codon_alleles($bvfoa);
        
        return 0 unless defined $ref_codon;
        
        return ( length($alt_codon) > length ($ref_codon) );
    }
    
    # structural variant
    elsif($bvfoa->base_variation_feature->isa('Bio::EnsEMBL::Variation::StructuralVariationFeature')) {
        
        # TO BE DONE, NO WAY OF KNOWING WHAT SEQUENCE IS INSERTED YET
        return 0;
        
        # must be an insertion
        return 0 unless insertion($bvfoa);
        
        my $bvfo = $bvfoa->base_variation_feature_overlap;
        
        my $cds_coords = $bvfo->cds_coords;
        
        if(scalar @$cds_coords == 1) {
            
            # wholly within exon
            if($cds_coords->[0]->isa('Bio::EnsEMBL::Mapper::Coordinate')) {
                1;
            }
        }
        
        # variant partially overlaps
        else {
            return 0;
        }
    }
    
    else {
        return 0;
    }
}

sub inframe_deletion {
    my $bvfoa = shift;
    
    # sequence variant
    if($bvfoa->base_variation_feature->isa('Bio::EnsEMBL::Variation::VariationFeature')) {
        my ($ref_codon, $alt_codon) = _get_codon_alleles($bvfoa);
        
        return 0 unless defined $ref_codon;
      
        return ( 
            (length($alt_codon) < length ($ref_codon)) &&
            ( ($ref_codon =~ /^\Q$alt_codon\E/) || ($ref_codon =~ /\Q$alt_codon\E$/) )
        );
    }
    
    # structural variant
    elsif($bvfoa->base_variation_feature->isa('Bio::EnsEMBL::Variation::StructuralVariationFeature')) {
        
        # must be a deletion
        return 0 unless deletion($bvfoa);
        
        my $bvfo       = $bvfoa->base_variation_feature_overlap;
        my $cds_coords = $bvfo->cds_coords;
        my $exons      = $bvfo->_exons;
        my $bvf        = $bvfo->base_variation_feature;
        
        # in exon
        return (
           scalar @$cds_coords == 1 and
           $cds_coords->[0]->isa('Bio::EnsEMBL::Mapper::Coordinate') and
           scalar grep {complete_within_feature($bvfoa, $_)} @$exons and
           $bvf->length() % 3 == 0
        );
    }
    
    else {
        return 0;
    }
}

sub stop_gained {
    my $bvfoa = shift;
    my $bvfo  = $bvfoa->base_variation_feature_overlap;
    
    return 0 unless $bvfoa->isa('Bio::EnsEMBL::Variation::TranscriptVariationAllele');

    my ($ref_pep, $alt_pep) = _get_peptide_alleles($bvfoa);
    
    return 0 unless defined $ref_pep;

    return ( ($alt_pep =~ /\*/) and ($ref_pep !~ /\*/) );
}

sub stop_lost {
    my $bvfoa = shift;
    
    # sequence variant
    if($bvfoa->isa('Bio::EnsEMBL::Variation::TranscriptVariationAllele')) {
        my ($ref_pep, $alt_pep) = _get_peptide_alleles($bvfoa);
        
        return 0 unless defined $ref_pep;
    
        return ( ($alt_pep !~ /\*/) and ($ref_pep =~ /\*/) );
    }
    
    # structural variant
    elsif($bvfoa->isa('Bio::EnsEMBL::Variation::TranscriptStructuralVariationAllele')) {
        return 0 unless deletion($bvfoa);
        
        my $tr  = $bvfoa->transcript;
        my $bvf = $bvfoa->base_variation_feature;
        
        my ($tr_crs, $tr_cre) = ($tr->coding_region_start, $tr->coding_region_end);
        return 0 unless defined($tr_crs) && defined($tr_cre);
        
        if($tr->strand == 1) {
            return overlap($tr_cre-2, $tr_cre, $bvf->start, $bvf->end);
        }
        else {
            return overlap($tr_crs, $tr_crs + 2, $bvf->start, $bvf->end);
        }
    }
    
    else {
        return 0;
    }
}

sub frameshift {
    my $bvfoa = shift;
    
    # sequence variant
    if($bvfoa->isa('Bio::EnsEMBL::Variation::TranscriptVariationAllele')) {

        return 0 if partial_codon($bvfoa);
    
        my $bvfo = $bvfoa->base_variation_feature_overlap;
    
        return 0 unless defined $bvfo->cds_start && defined $bvfo->cds_end;
        
        my $var_len = $bvfo->cds_end - $bvfo->cds_start + 1;
    
        my $allele_len = $bvfoa->seq_length;
    
        # if the allele length is undefined then we can't call a frameshift
    
        return 0 unless defined $allele_len;
    
        return abs( $allele_len - $var_len ) % 3;
    }
    
    # structural variant
    elsif($bvfoa->isa('Bio::EnsEMBL::Variation::TranscriptStructuralVariationAllele')) {
        my $bvf   = $bvfoa->base_variation_feature;
        my $exons = $bvfoa->base_variation_feature_overlap->_exons;
        
        return (
            (
                deletion($bvfoa) or
                copy_number_loss($bvfoa)
            ) and
            scalar grep {complete_within_feature($bvfoa, $_)} @$exons and
            $bvf->length % 3 != 0
        );
        
        # TODO INSERTIONS
    }
    
    else {
        return 0;
    }
}

sub partial_codon {
    my $bvfoa = shift;
    
    my $bvfo = $bvfoa->base_variation_feature_overlap;
    
    return 0 unless defined $bvfo->translation_start;

    my $cds_length = length $bvfo->_translateable_seq;

    my $codon_cds_start = ($bvfo->translation_start * 3) - 2;

    my $last_codon_length = $cds_length - ($codon_cds_start - 1);
    
    return ( $last_codon_length < 3 and $last_codon_length > 0 );
}

sub coding_unknown {
    my $bvfoa = shift;
    
    # sequence variant
    if($bvfoa->isa('Bio::EnsEMBL::Variation::TranscriptVariationAllele')) {
        return (within_cds($bvfoa) and ((not $bvfoa->peptide) or (not _get_peptide_alleles($bvfoa)) or ($bvfoa->peptide =~ /X/)) and (not (frameshift($bvfoa) or inframe_deletion($bvfoa))));
    }
    
    # structural variant
    elsif($bvfoa->isa('Bio::EnsEMBL::Variation::TranscriptStructuralVariationAllele')) {
        return (within_cds($bvfoa) and not(inframe_insertion($bvfoa) or inframe_deletion($bvfoa) or frameshift($bvfoa)));
    }
    
    else {
        return 0;
    }
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

    my $bvf  = $vfo->base_variation_feature;
    my $feat = $vfo->feature;

    return ( ($bvf->start <= $feat->start) && ($bvf->end >= $feat->end) ); 
}

1;

