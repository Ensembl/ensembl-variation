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

=head1 NAME

Bio::EnsEMBL::Variation::TranscriptStructuralVariationAllele

=head1 SYNOPSIS

    use Bio::EnsEMBL::Variation::TranscriptStructuralVariationAllele;

=head1 DESCRIPTION

A TranscriptStructuralVariationAllele object represents a single allele of a
TranscriptStructuralVariation.

=cut

package Bio::EnsEMBL::Variation::TranscriptStructuralVariationAllele;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::StructuralVariationOverlapAllele Bio::EnsEMBL::Variation::BaseTranscriptVariationAllele);

sub new_fast {
    my ($class, $hashref) = @_;
    
    # swap a transcript_structural_variation argument for a structural_variation_overlap one

    if ($hashref->{transcript_structural_variation}) {
        $hashref->{structural_variation_overlap} = delete $hashref->{transcript_structural_variation};
    }
    
    # and call the superclass

    return $class->SUPER::new_fast($hashref);
}

=head2 transcript_structural_variation

  Description: Get the associated TranscriptStructuralVariation
  Returntype : Bio::EnsEMBL::Variation::TranscriptStructuralVariation
  Exceptions : none
  Status     : Stable

=cut

sub transcript_structural_variation {
    my ($self, $svo) = @_;
    if ($svo) {
        assert_ref($svo, 'Bio::EnsEMBL::Variation::TranscriptStructuralVariation');
    }
    return $self->base_variation_feature_overlap($svo);
}


1;

