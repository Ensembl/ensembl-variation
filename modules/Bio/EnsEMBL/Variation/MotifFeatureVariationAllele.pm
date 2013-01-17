=head1 LICENSE

 Copyright (c) 1999-2013 The European Bioinformatics Institute and
 Genome Research Limited.  All rights reserved.

 This software is distributed under a modified Apache license.
 For license details, please see

   http://www.ensembl.org/info/about/legal/code_licence.html

=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <dev@ensembl.org>.

 Questions may also be sent to the Ensembl help desk at
 <helpdesk@ensembl.org>.

=cut

package Bio::EnsEMBL::Variation::MotifFeatureVariationAllele;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(overlap);
use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);

use base qw(Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele);

sub new_fast {
    my ($self, $hashref) = @_;
    
    # swap a motif_feature_variation argument for a variation_feature_overlap one

    if ($hashref->{motif_feature_variation}) {
        $hashref->{variation_feature_overlap} = delete $hashref->{motif_feature_variation};
    }
    
    # and call the superclass

    return $self->SUPER::new_fast($hashref);
}

=head2 motif_feature_variation

  Description: Get/set the associated MotifFeatureVariation
  Returntype : Bio::EnsEMBL::Variation::MotifFeatureVariation
  Status     : Stable

=cut

sub motif_feature_variation {
    my ($self, $mfv) = shift;
    assert_ref($mfv, 'Bio::EnsEMBL::Variation::MotifFeatureVariation') if $mfv;
    return $self->variation_feature_overlap($mfv);
}

=head2 motif_feature

  Description: Get/set the associated MotifFeature
  Returntype : Bio::EnsEMBL::Funcgen::MotifFeature
  Status     : Stable

=cut

sub motif_feature {
    my $self = shift;
    return $self->motif_feature_variation->motif_feature;
}

=head2 motif_start

  Description: Get the (1-based) relative position of the variation feature start in the motif
  Returntype : integer
  Status     : At Risk

=cut

sub motif_start {

    my $self = shift;
    unless ($self->{motif_start}) {
        my $mf = $self->motif_feature;
        my $vf = $self->variation_feature;
        
        return undef unless defined $vf->seq_region_start && defined $mf->seq_region_start;
       
        my $mf_start = $vf->seq_region_start - $mf->seq_region_start + 1;

        # adjust if the motif is on the reverse strand

        $mf_start = $mf->binding_matrix->length - $mf_start + 1 if $mf->strand < 0;
        
        # check that we're in bounds

        return undef if $mf_start > $mf->length;
        $self->{motif_start} = $mf_start;
    }
    return $self->{motif_start};
}

=head2 motif_end

  Description: Get the (1-based) relative position of the variation feature end in the motif
  Returntype : integer
  Status     : At Risk

=cut

sub motif_end {

    my $self = shift;
    unless ($self->{motif_end}) {
        my $mf = $self->motif_feature;
        my $vf = $self->variation_feature;
        
        return undef unless defined $vf->seq_region_end && defined $mf->seq_region_start;
       
        my $mf_end = $vf->seq_region_end - $mf->seq_region_start + 1;

        # adjust if the motif is on the reverse strand

        $mf_end = $mf->binding_matrix->length - $mf_end + 1 if $mf->strand < 0;
        
        # check that we're in bounds

        return undef if $mf_end < 1;
        $self->{motif_end} = $mf_end;
    }
    return $self->{motif_end};
}

=head2 in_informative_position

  Description: Identify if the variation feature falls in a high information position of the motif
  Returntype : boolean
  Status     : At Risk

=cut

sub in_informative_position {
    my $self = shift;

    # we can only call this for true SNPs
    unless ($self->{in_informative_position}) {
        my $vf = $self->variation_feature;

        unless (($vf->start == $vf->end) && ($self->variation_feature_seq ne '-')) {
            return undef;
        }

        # get the 1-based position

        my $start = $self->motif_start;

        return undef unless defined $start && $start >= 1 && $start <= $self->motif_feature->length;

        $self->{in_informative_position} = $self->motif_feature->binding_matrix->is_position_informative($start);
    }
    return $self->{in_informative_position};
}

=head2 motif_score_delta

  Description: Calculate the difference in motif score between the reference and variant sequences
  Returntype : float
  Status     : At Risk

=cut

sub motif_score_delta {

    my $self    = shift;
    my $linear  = shift;
    
    unless ($self->{motif_score_delta}) {

        my $vf = $self->motif_feature_variation->variation_feature;
        my $mf = $self->motif_feature;

        my $allele_seq      = $self->feature_seq;
        my $ref_allele_seq  = $self->motif_feature_variation->get_reference_MotifFeatureVariationAllele->feature_seq;

        if ($allele_seq eq '-' || 
            $ref_allele_seq eq '-' || 
            length($allele_seq) != length($ref_allele_seq)) {
            # we can't call a score because the sequence will change length
            return undef;
        }

        my ($mf_start, $mf_end) = ($self->motif_start, $self->motif_end);

        return undef unless defined $mf_start && defined $mf_end;

        my $mf_seq = $self->motif_feature_variation->_motif_feature_seq;
        my $mf_seq_length = length($mf_seq);
        
        # trim allele_seq
        if($mf_start < 1) {
            $allele_seq = substr($allele_seq, 1 - $mf_start);
            $mf_start = 1;
        }
        
        if($mf_end > $mf_seq_length) {
            $allele_seq = substr($allele_seq, 0, $mf_seq_length - $mf_start + 1);
            $mf_end = $mf_seq_length;
        }

        my $var_len = length($allele_seq);

        return undef if $var_len > $mf->length;

        my $matrix = $mf->binding_matrix;

        # get the binding affinity of the reference sequence
        my $ref_affinity = $matrix->relative_affinity($mf_seq, $linear);

        # splice in the variant sequence (0-based)
        substr($mf_seq, $mf_start - 1, $var_len) = $allele_seq;
        
        # check length hasn't changed
        return undef if length($mf_seq) != $mf_seq_length;

        # and get the affinity of the variant sequence
        my $var_affinity = $matrix->relative_affinity($mf_seq, $linear);

        $self->{motif_score_delta} = ($var_affinity - $ref_affinity);
    }

    return $self->{motif_score_delta};
}

1;
