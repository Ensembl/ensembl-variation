=head1 LICENSE

 Copyright (c) 1999-2013 The European Bioinformatics Institute and
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

package Bio::EnsEMBL::Variation::StructuralVariationOverlapAllele;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::BaseVariationFeatureOverlapAllele);

sub new_fast {
    my ($class, $hashref) = @_;

    # swap a transcript_variation argument for a variation_feature_overlap one

    if ($hashref->{structural_variation_overlap}) {
        $hashref->{base_variation_feature_overlap} = 
            delete $hashref->{structural_variation_overlap};
    }
    
    # and call the superclass

    return $class->SUPER::new_fast($hashref);
}

=head2 structural_variation_overlap

  Description: Get the associated StructuralVariationOverlap
  Returntype : Bio::EnsEMBL::Variation::StructuralVariationOverlap
  Exceptions : none
  Status     : Stable

=cut

sub structural_variation_overlap {
    my ($self, $svo) = @_;
    if ($svo) {
        assert_ref($svo, 'Bio::EnsEMBL::Variation::StructuralVariationOverlap');
    }
    return $self->base_variation_feature_overlap($svo);
}


=head2 structural_variation_feature

  Description: Get the associated StructuralVariationFeature
  Returntype : Bio::EnsEMBL::Variation::StructuralVariationFeature
  Exceptions : none
  Status     : Stable

=cut

sub structural_variation_feature {
    my $self = shift;
    return $self->structural_variation_overlap->structural_variation_feature;
}

1;

