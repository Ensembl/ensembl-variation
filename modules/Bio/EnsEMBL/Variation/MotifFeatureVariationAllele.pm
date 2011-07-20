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

package Bio::EnsEMBL::Variation::MotifFeatureVariationAllele;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(overlap);

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

sub motif_feature_variation {
    my $self = shift;
    return $self->variation_feature_overlap(@_);
}

sub motif_feature {
    my $self = shift;
    return $self->motif_feature_variation->motif_feature;
}

sub in_informative_position {
    my $self = shift;

    my $mf = $self->motif_feature;
    my $vf = $self->variation_feature;

    # we can only call this for true SNPs

    unless (($vf->start == $vf->end) && ($self->variation_feature_seq ne '-')) {
        return undef;
    }
    
    # get the (1-based) relative start of the vf with respect to this mf
    
    my $mf_start = $vf->seq_region_start - $mf->seq_region_start + 1;

    # check that we're in bounds

    return undef if ( ($mf_start < 1) || ($mf_start > $mf->length) );
        
    return $mf->is_position_informative($mf_start);
}

sub binding_affinity_change {
    my $self    = shift;
    my $linear  = shift;
    
    unless ($self->{binding_affinity_change}) {
        
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
        
        # get the relative start of the vf with respect to this mf
    
        my $mf_start = $vf->seq_region_start - $mf->seq_region_start;
        
        return undef if $mf_start < 0;
        
        my $var_len = length($allele_seq);
        
        return undef if $var_len > $mf->length;
        
        my $mf_seq = $mf->seq;
        
        my $matrix = $mf->binding_matrix;
        
        # get the binding affinity of the reference sequence
        my $ref_affinity = $matrix->relative_affinity($mf_seq, $linear);
        
        # splice in the variant sequence
        substr($mf_seq, $mf_start, $var_len) = $allele_seq;
        
        # and get the affinity of the variant sequence
        my $var_affinity = $matrix->relative_affinity($mf_seq, $linear);
        
        $self->{binding_affinity_change} = ($var_affinity - $ref_affinity);
    }
    
    return $self->{binding_affinity_change};
}

1;
