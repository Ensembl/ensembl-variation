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

sub motif_feature_variation {
    my $self = shift;
    return $self->variation_feature_overlap;
}

sub motif_feature {
    my $self = shift;
    return $self->variation_feature_overlap->feature;
}

sub binding_affinity_change {
    my $self    = shift;
    my $linear  = shift;
    
    unless ($self->{binding_affinity_change}) {
        
        my $vf = $self->motif_feature_variation->variation_feature;
        my $mf = $self->motif_feature_variation->motif_feature;;
        
        my $allele_seq      = $self->feature_seq;
        my $ref_allele_seq  = $self->motif_feature_variation->reference_allele->feature_seq;
        
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
        
        warn "seq: $mf_seq\n";
        
        # splice in the variant sequence
        substr($mf_seq, $mf_start, $var_len) = $allele_seq;
        
        warn "seq: $mf_seq\n";
        
        # and get the affinity of the variant sequence
        my $var_affinity = $matrix->relative_affinity($mf_seq, $linear);
        
        warn "ref aff: $ref_affinity var aff: $var_affinity\n";
        
        $self->{binding_affinity_change} = ($var_affinity - $ref_affinity);
    }
    
    return $self->{binding_affinity_change};
}

1;
