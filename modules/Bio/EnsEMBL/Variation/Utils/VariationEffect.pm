=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file excepst in compliance with the License.
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
 developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

 Questions may also be sent to the Ensembl help desk at
 <http://www.ensembl.org/Help/Contact>.

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

our @EXPORT_OK = qw(overlap within_feature within_cds MAX_DISTANCE_FROM_TRANSCRIPT within_intron stop_lost stop_retained affects_start_codon frameshift $UPSTREAM_DISTANCE $DOWNSTREAM_DISTANCE);

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
    my ($bvfoa, $feat, $bvfo, $bvf) = @_;
    $bvf ||= $bvfoa->base_variation_feature;
    
    # sequence variant will have alleles
    if($bvf->isa('Bio::EnsEMBL::Variation::VariationFeature')) {
        my ($ref_allele, $alt_allele) = _get_alleles(@_);
        return (
            (defined($ref_allele) && ($alt_allele eq '' || length($alt_allele) < length($ref_allele)) and $ref_allele) or
            $bvf->allele_string =~ /deletion/i
        );
    }
    
    # structural variant depends on class
    if($bvf->isa('Bio::EnsEMBL::Variation::StructuralVariationFeature')) {
        return (
            ($bvf->class_SO_term(undef, 1) eq 'deletion') or
            ($bvf->class_SO_term(undef, 1) =~ /deletion/i) or
            ($bvf->class_SO_term(undef, 1) =~ /loss/i)
        );
    }
    
    else { return 0; }
}

sub insertion {
    my ($bvfoa, $feat, $bvfo, $bvf) = @_;
    $bvf ||= $bvfoa->base_variation_feature;
    
    # sequence variant will have alleles
    if($bvf->isa('Bio::EnsEMBL::Variation::VariationFeature')) {
        my ($ref_allele, $alt_allele) = _get_alleles(@_);
        return (
            (defined($ref_allele) && ($ref_allele eq '' || length($alt_allele) > length($ref_allele)) and $alt_allele) or
            $bvf->allele_string =~ /insertion/i
        );
    }
    
    # structural variant depends on class
    if($bvf->isa('Bio::EnsEMBL::Variation::StructuralVariationFeature')) {
        return (
            duplication(@_) or
            tandem_duplication(@_) or
            ($bvf->class_SO_term(undef, 1) eq 'insertion') or
            ($bvf->class_SO_term(undef, 1) =~ /insertion/i) or
            ($bvf->class_SO_term(undef, 1) =~ /gain/i)
        );
    }
    
    else { return 0; }
}

sub copy_number_gain {
    my ($bvfoa, $feat, $bvfo, $bvf) = @_;
    $bvf ||= $bvfoa->base_variation_feature;
    
    return (duplication(@_) or tandem_duplication(@_) or $bvf->class_SO_term(undef, 1) =~ /gain/i);
}

sub copy_number_loss {
    my ($bvfoa, $feat, $bvfo, $bvf) = @_;
    $bvf ||= $bvfoa->base_variation_feature;
    
    return $bvf->class_SO_term(undef, 1) =~ /loss/i;
}

sub duplication {
    my ($bvfoa, $feat, $bvfo, $bvf) = @_;
    $bvf ||= $bvfoa->base_variation_feature;
    
    return (
        (
            ($bvf->class_SO_term(undef, 1) eq 'duplication') or
            ($bvf->class_SO_term(undef, 1) =~ /duplication/i)
        ) and
        (not tandem_duplication(@_))
    );
}

sub tandem_duplication {
    my ($bvfoa, $feat, $bvfo, $bvf) = @_;
    $bvf ||= $bvfoa->base_variation_feature;
    
    # for sequence variants, check sequence vs ref
    if($bvf->isa('Bio::EnsEMBL::Variation::VariationFeature')) {
        my ($ref_allele, $alt_allele) = _get_alleles(@_);
        
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
            ($bvf->class_SO_term(undef, 1) eq 'tandem_duplication') or
            ($bvf->class_SO_term(undef, 1) =~ /tandem_duplication/i) 
        );
    }
}

sub feature_ablation {
    my ($bvfoa, $feat, $bvfo, $bvf) = @_;
    $feat ||= $bvfoa->base_variation_feature_overlap->feature;
    
    return (complete_overlap_feature($bvfoa, $feat) and deletion(@_));
}

sub feature_amplification {
    my ($bvfoa, $feat, $bvfo, $bvf) = @_;
    $feat ||= $bvfoa->base_variation_feature_overlap->feature;
    
    return (complete_overlap_feature($bvfoa, $feat) and copy_number_gain(@_));
}

sub feature_elongation {
    my ($bvfoa, $feat, $bvfo, $bvf) = @_;
    $feat ||= $bvfoa->base_variation_feature_overlap->feature;
    
    return 0 if $bvfoa->isa('Bio::EnsEMBL::Variation::TranscriptVariationAllele');
    
    return (
        complete_within_feature($bvfoa, $feat) and
        (copy_number_gain(@_) or insertion(@_))
    );
}

sub feature_truncation {
    my ($bvfoa, $feat, $bvfo, $bvf) = @_;
    $feat ||= $bvfoa->base_variation_feature_overlap->feature;
    
    return 0 if $bvfoa->isa('Bio::EnsEMBL::Variation::TranscriptVariationAllele');
    
    return (
        (partial_overlap_feature($bvfoa, $feat) or complete_within_feature($bvfoa, $feat)) and
        (copy_number_loss(@_) or deletion(@_))
    );
}

#sub transcript_fusion {
#    #my ($bvfoa, $feat, $bvfo, $bvf) = @_;
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
    my ($bvfoa, $feat, $bvfo, $bvf) = @_;
    $bvf  ||= $bvfoa->base_variation_feature;
    $feat ||= $bvfoa->base_variation_feature_overlap->feature;

    return _upstream($bvf, $feat, $UPSTREAM_DISTANCE);
}

sub downstream {
    my ($bvfoa, $feat, $bvfo, $bvf) = @_;
    $bvf  ||= $bvfoa->base_variation_feature;
    $feat ||= $bvfoa->base_variation_feature_overlap->feature;

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
    return within_feature(@_);
}

sub within_nmd_transcript {
    my ($bvfoa, $feat, $bvfo, $bvf) = @_;
    $feat ||= $bvfoa->base_variation_feature_overlap->feature;

    return ( within_transcript(@_) and ($feat->biotype eq 'nonsense_mediated_decay') );
}

sub within_non_coding_gene {
    my ($bvfoa, $feat, $bvfo, $bvf) = @_;
    $feat ||= $bvfoa->base_variation_feature_overlap->feature;
    
    return ( within_transcript(@_) and (not $feat->translation) and (not within_mature_miRNA(@_)));
}

sub non_coding_exon_variant {
    my ($bvfoa, $feat, $bvfo, $bvf) = @_;
    $bvfo ||= $bvfoa->base_variation_feature_overlap;
    
    return 0 unless within_non_coding_gene(@_);
    
    my $exons = $bvfo->_exons;
    
    if(scalar grep {overlap($bvf->start, $bvf->end, $_->start, $_->end)} @$exons) {
        return 1;
    }
    else {
        return 0;
    }
}

sub within_miRNA {
    my ($bvfoa, $feat, $bvfo, $bvf) = @_;
    
    # don't call this for now
    
    return 0;
    $feat ||= $bvfoa->base_variation_feature_overlap->feature;
    
    return ( within_transcript(@_) and ($feat->biotype eq 'miRNA') );
}

sub within_mature_miRNA {
    my ($bvfoa, $feat, $bvfo, $bvf) = @_;
    $bvfo ||= $bvfoa->base_variation_feature_overlap;
    $bvf  ||= $bvfo->base_variation_feature;
    $feat ||= $bvfo->feature;
        
    return 0 unless ( within_transcript(@_) and ($feat->biotype eq 'miRNA') );
        
    my ($attribute) = @{ $feat->get_all_Attributes('miRNA') };
    
    if (defined $attribute && $attribute->value =~ /(\d+)-(\d+)/) { 
        for my $coord ($bvfo->_mapper->cdna2genomic($1, $2, $feat->strand)) {
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
    my ($bvfoa, $feat, $bvfo, $bvf) = @_;
    $bvfo ||= $bvfoa->base_variation_feature_overlap;
    $feat ||= $bvfo->feature;
    
    return $feat->strand == 1 ? 
        $bvfo->_intron_effects->{start_splice_site} :
        $bvfo->_intron_effects->{end_splice_site};
}

sub acceptor_splice_site {
    my ($bvfoa, $feat, $bvfo, $bvf) = @_;
    $bvfo ||= $bvfoa->base_variation_feature_overlap;
    $feat ||= $bvfo->feature;
    
    return $feat->strand == 1 ? 
        $bvfo->_intron_effects->{end_splice_site} :
        $bvfo->_intron_effects->{start_splice_site};
}

sub essential_splice_site {
    return ( acceptor_splice_site(@_) or donor_splice_site(@_) );
}

sub splice_region {
    my ($bvfoa, $feat, $bvfo, $bvf) = @_;
    $bvfo ||= $bvfoa->base_variation_feature_overlap;
    
    return 0 if donor_splice_site(@_);
    return 0 if acceptor_splice_site(@_);
    return 0 if essential_splice_site(@_);

    return $bvfo->_intron_effects->{splice_region};
}

sub within_intron {
    my ($bvfoa, $feat, $bvfo, $bvf) = @_;
    $bvfo ||= $bvfoa->base_variation_feature_overlap;

    return $bvfo->_intron_effects->{intronic};
}

sub within_cds {
    my ($bvfoa, $feat, $bvfo, $bvf) = @_;
    $bvfo ||= $bvfoa->base_variation_feature_overlap;
    $feat ||= $bvfo->feature;
    $bvf  ||= $bvfo->base_variation_feature;
    
    my $cds_coords = $bvfo->cds_coords;
    
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

    if (defined $feat->translation &&
        $bvfo->_intron_effects->{within_frameshift_intron}) {
 
        return overlap(
            $bvf->start, 
            $bvf->end, 
            $feat->coding_region_start,
            $feat->coding_region_end,
        );
    }
        
    return 0;
}

sub within_cdna {
    my ($bvfoa, $feat, $bvfo, $bvf) = @_;
    $bvfo ||= $bvfoa->base_variation_feature_overlap;
    $feat ||= $bvfo->feature;
    
    my $cdna_coords = $bvfo->cdna_coords;
    
    if (@$cdna_coords > 0) {
        for my $coord (@$cdna_coords) {
            if ($coord->isa('Bio::EnsEMBL::Mapper::Coordinate')) {
                if ($coord->end > 0 && $coord->start <= $feat->length) {
                    return 1;
                }
            }
        }
    }
    
    # we also need to check if the vf is in a frameshift intron within the cDNA

    if ($bvfo->_intron_effects->{within_frameshift_intron}) {
        return within_transcript(@_); 
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
    my ($bvfoa, $feat, $bvfo, $bvf) = @_;
    $bvfo ||= $bvfoa->base_variation_feature_overlap;
    $feat ||= $bvfo->feature;
    $bvf  ||= $bvfo->base_variation_feature;
    
    my $five_prime_of_coding = 
        $feat->strand == 1 ? 
            _before_coding($bvf, $feat) : 
            _after_coding($bvf, $feat);
    
    return ( $five_prime_of_coding and within_cdna(@_) );
}

sub within_3_prime_utr {
    my ($bvfoa, $feat, $bvfo, $bvf) = @_;
    $bvfo ||= $bvfoa->base_variation_feature_overlap;
    $feat ||= $bvfo->feature;
    $bvf  ||= $bvfo->base_variation_feature;
    
    my $three_prime_of_coding = 
        $feat->strand == 1 ? 
            _after_coding($bvf, $feat) : 
            _before_coding($bvf, $feat);
    
    return ( $three_prime_of_coding and within_cdna(@_) );
}

sub complex_indel {
    my ($bvfoa, $feat, $bvfo, $bvf) = @_;
    $bvfo ||= $bvfoa->base_variation_feature_overlap;
    $bvf  ||= $bvfo->base_variation_feature;
   
    # pass the no_db flag to var_class to ensure we don't rely on the database for it 
    # as it may not have been set at this stage in the pipeline
    my $class = $bvf->var_class(1);

    return 0 unless $class =~ /insertion|deletion|indel/;

    return @{ $bvfo->cds_coords } > 1;
}

sub _get_peptide_alleles {
    my ($bvfoa, $feat, $bvfo, $bvf) = @_;

    return () unless defined $bvfoa;
    #return () if frameshift(@_);

    my $alt_pep = $bvfoa->peptide;
    
    return () unless defined $alt_pep;
    
    my $ref_pep = _get_ref_pep(@_);
    
    return () unless defined $ref_pep;
    
    $ref_pep = '' if $ref_pep eq '-';
    $alt_pep = '' if $alt_pep eq '-';
    
    return ($ref_pep, $alt_pep);
}

sub _get_ref_pep {
    my ($bvfoa, $feat, $bvfo, $bvf) = @_;
    $bvfo ||= $bvfoa->base_variation_feature_overlap;
    return $bvfo->get_reference_TranscriptVariationAllele->peptide;
}

sub _get_codon_alleles {
    my ($bvfoa, $feat, $bvfo, $bvf) = @_;
    
    return () if frameshift(@_);

    my $alt_codon = $bvfoa->codon;
    
    return () unless defined $alt_codon;
    
    my $ref_codon = $bvfo->get_reference_TranscriptVariationAllele->codon;
    
    return () unless defined $ref_codon;
    
    $ref_codon = '' if $ref_codon eq '-';
    $alt_codon = '' if $alt_codon eq '-';
    
    return ($ref_codon, $alt_codon);
}

sub _get_alleles {
    my ($bvfoa, $feat, $bvfo, $bvf) = @_;
    $bvfo ||= $bvfoa->base_variation_feature_overlap;
    
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
    return 0 unless defined $alt_pep && $alt_pep =~/^\*/; 

    ## handle inframe insertion of a stop just before the stop (no ref peptide)
    if( $bvfoa->isa('Bio::EnsEMBL::Variation::TranscriptVariationAllele') &&
        defined $bvfoa->transcript_variation->transcript->translation &&
         $bvfoa->transcript_variation->translation_start() > 
           $bvfoa->transcript_variation->transcript->translation->length()
	){

	return 1;
    }
    return 0 unless $ref_pep;

    return ( $alt_pep =~ /^\*/ && $ref_pep =~ /^\*/ );
}

sub affects_start_codon {
    my ($bvfoa, $feat, $bvfo, $bvf) = @_;
    $bvfo ||= $bvfoa->base_variation_feature_overlap;
    $feat ||= $bvfo->feature;
    $bvf  ||= $bvfo->base_variation_feature;
    
    # sequence variant
    if($bvfo->isa('Bio::EnsEMBL::Variation::TranscriptVariation')) {
        my ($ref_pep, $alt_pep) = _get_peptide_alleles(@_);
    
        return 0 unless $ref_pep;
    
        return ( ($bvfo->translation_start == 1) and (substr($ref_pep,0,1) ne substr($alt_pep,0,1)) );
    }
    
    # structural variant
    elsif($bvfo->isa('Bio::EnsEMBL::Variation::TranscriptStructuralVariation')) {        
        my ($tr_crs, $tr_cre) = ($feat->coding_region_start, $feat->coding_region_end);
        return 0 unless defined($tr_crs) && defined($tr_cre);
        
        if($feat->strand == 1) {
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
    my ($ref_pep, $alt_pep) = _get_peptide_alleles(@_);
    
    return 0 unless $ref_pep;

    return ( ($alt_pep eq $ref_pep) and (not stop_retained(@_) and ($alt_pep !~ /X/) and ($ref_pep !~ /X/)) );
}

sub missense_variant {
    my ($ref_pep, $alt_pep) = _get_peptide_alleles(@_);
    
    return 0 unless defined $ref_pep;
    
    return 0 if affects_start_codon(@_);
    return 0 if stop_lost(@_);
    return 0 if stop_gained(@_);
    return 0 if partial_codon(@_);
    return 0 if frameshift(@_);

    return 0 if inframe_deletion(@_);
    return 0 if inframe_insertion(@_);
    
    return ( $ref_pep ne $alt_pep ) && ( length($ref_pep) == length($alt_pep) );
}

sub inframe_insertion {
    my ($bvfoa, $feat, $bvfo, $bvf) = @_;
    $bvfo ||= $bvfoa->base_variation_feature_overlap;
    $bvf  ||= $bvfo->base_variation_feature;

    # sequence variant
    if($bvf->isa('Bio::EnsEMBL::Variation::VariationFeature')) {
        my ($ref_codon, $alt_codon) = _get_codon_alleles(@_);

        return 0 if affects_start_codon(@_);
        return 0 unless defined $ref_codon;

        return ( length($alt_codon) > length ($ref_codon) );
    }
    
    # structural variant
    elsif($bvf->isa('Bio::EnsEMBL::Variation::StructuralVariationFeature')) {

        # TO BE DONE, NO WAY OF KNOWING WHAT SEQUENCE IS INSERTED YET
        return 0;
        
        # must be an insertion
        return 0 unless insertion(@_);
        
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
    my ($bvfoa, $feat, $bvfo, $bvf) = @_;
    $bvfo ||= $bvfoa->base_variation_feature_overlap;
    $bvf  ||= $bvfo->base_variation_feature;
    
    # sequence variant
    if($bvf->isa('Bio::EnsEMBL::Variation::VariationFeature')) {
        my ($ref_codon, $alt_codon) = _get_codon_alleles(@_);
        
        return 0 unless defined $ref_codon;
        return 0 unless length($alt_codon) < length ($ref_codon);
        
        # simple string match
        return 1 if ($ref_codon =~ /^\Q$alt_codon\E/) || ($ref_codon =~ /\Q$alt_codon\E$/);
        
        # try a more complex string match; matching part may be in the middle
        # first trim matching bases from start of string
        while($ref_codon && $alt_codon && substr($ref_codon, 0, 1) eq substr($alt_codon, 0, 1)) {
          $ref_codon = substr($ref_codon, 1);
          $alt_codon = substr($alt_codon, 1);
        }
        
        # now trim ends
        while($ref_codon && $alt_codon && substr($ref_codon, -1, 1) eq substr($alt_codon, -1, 1)) {
          $ref_codon = substr($ref_codon, 0, length($ref_codon) - 1);
          $alt_codon = substr($alt_codon, 0, length($alt_codon) - 1);
        }
        
        # if nothing remains of $alt_codon,
        # then it fully matched a part in the middle of $ref_codon
        return length($alt_codon) == 0;
    }
    
    # structural variant
    elsif($bvf->isa('Bio::EnsEMBL::Variation::StructuralVariationFeature')) {
        
        # must be a deletion
        return 0 unless deletion(@_);
        
        my $cds_coords = $bvfo->cds_coords;
        my $exons      = $bvfo->_exons;
        
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
    my ($bvfoa, $feat, $bvfo, $bvf) = @_;
    
    return 0 unless $bvfoa->isa('Bio::EnsEMBL::Variation::TranscriptVariationAllele');
    ## check for inframe insertion before stop 
    return 0 if stop_retained($bvfoa);

    my ($ref_pep, $alt_pep) = _get_peptide_alleles(@_);
    
    return 0 unless defined $ref_pep;

    return ( ($alt_pep =~ /\*/) and ($ref_pep !~ /\*/) );
}

sub stop_lost {
    my ($bvfoa, $feat, $bvfo, $bvf) = @_;
    $bvfo ||= $bvfoa->base_variation_feature_overlap;
    $bvf  ||= $bvfo->base_variation_feature;
    $feat ||= $bvfo->feature;
    
    # sequence variant
    if($bvfoa->isa('Bio::EnsEMBL::Variation::TranscriptVariationAllele')) {
        
        # special case frameshift
#        if(frameshift(@_)) {
#          my $ref_pep = _get_ref_pep(@_);
#          return $ref_pep && $ref_pep =~ /\*/;
#        }
        
        my ($ref_pep, $alt_pep) = _get_peptide_alleles(@_);
        
        return 0 unless defined $ref_pep;
    
        return ( ($alt_pep !~ /\*/) and ($ref_pep =~ /\*/) );
    }
    
    # structural variant
    elsif($bvfoa->isa('Bio::EnsEMBL::Variation::TranscriptStructuralVariationAllele')) {
        return 0 unless deletion(@_);
        
        my ($tr_crs, $tr_cre) = ($feat->coding_region_start, $feat->coding_region_end);
        return 0 unless defined($tr_crs) && defined($tr_cre);
        
        if($feat->strand == 1) {
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
    my ($bvfoa, $feat, $bvfo, $bvf) = @_;
    $bvfo ||= $bvfoa->base_variation_feature_overlap;
    
    # sequence variant
    if($bvfoa->isa('Bio::EnsEMBL::Variation::TranscriptVariationAllele')) {

        return 0 if partial_codon(@_);
    
        return 0 unless defined $bvfo->cds_start && defined $bvfo->cds_end;
        
        my $var_len = $bvfo->cds_end - $bvfo->cds_start + 1;
    
        my $allele_len = $bvfoa->seq_length;
    
        # if the allele length is undefined then we can't call a frameshift
    
        return 0 unless defined $allele_len;
    
        return abs( $allele_len - $var_len ) % 3;
    }
    
    # structural variant
    elsif($bvfoa->isa('Bio::EnsEMBL::Variation::TranscriptStructuralVariationAllele')) {
        my $exons = $bvfo->_exons;
        
        return (
            (
                deletion(@_) or
                copy_number_loss(@_)
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
    my ($bvfoa, $feat, $bvfo, $bvf) = @_;
    $bvfo ||= $bvfoa->base_variation_feature_overlap;
    
    return 0 unless defined $bvfo->translation_start;

    my $cds_length = length $bvfo->_translateable_seq;

    my $codon_cds_start = ($bvfo->translation_start * 3) - 2;

    my $last_codon_length = $cds_length - ($codon_cds_start - 1);
    
    return ( $last_codon_length < 3 and $last_codon_length > 0 );
}

sub coding_unknown {
    my ($bvfoa, $feat, $bvfo, $bvf) = @_;
    
    # sequence variant
    if($bvfoa->isa('Bio::EnsEMBL::Variation::TranscriptVariationAllele')) {
        return (within_cds(@_) and ((not $bvfoa->peptide) or (not _get_peptide_alleles(@_)) or ($bvfoa->peptide =~ /X/)) and (not (frameshift(@_) or inframe_deletion(@_))));
    }
    
    # structural variant
    elsif($bvfoa->isa('Bio::EnsEMBL::Variation::TranscriptStructuralVariationAllele')) {
        return (within_cds(@_) and not(inframe_insertion(@_) or inframe_deletion(@_) or frameshift(@_)));
    }
    
    else {
        return 0;
    }
}

#package Bio::EnsEMBL::Variation::RegulatoryFeatureVariationAllele;

sub within_regulatory_feature {
    my ($bvfoa, $feat, $bvfo, $bvf) = @_;
    return within_feature($bvfoa, $feat);
}

#package Bio::EnsEMBL::Variation::ExternalFeatureVariationAllele;

sub within_external_feature {
    my ($bvfoa, $feat, $bvfo, $bvf) = @_;
    return (within_feature($bvfoa, $feat) and (not within_miRNA_target_site(@_)));
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
    my ($bvfoa, $feat, $bvfo, $bvf) = @_;
    return within_feature($bvfoa, $feat);
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

