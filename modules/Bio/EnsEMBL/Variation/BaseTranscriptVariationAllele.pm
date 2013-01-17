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

Bio::EnsEMBL::Variation::BaseTranscriptVariationAllele

=head1 SYNOPSIS

    use Bio::EnsEMBL::Variation::BaseTranscriptVariationAllele;

=head1 DESCRIPTION

An helper class for representing an overlap of a Transcript and a
Variation allele (either sequence or structural). Should not be invoked
directly.

=cut

package Bio::EnsEMBL::Variation::BaseTranscriptVariationAllele;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele);

=head2 base_transcript_variation

  Description: Get/set the associated BaseTranscriptVariation
  Returntype : Bio::EnsEMBL::Variation::BaseTranscriptVariation
  Exceptions : throws if the argument is the wrong type
  Status     : Stable

=cut

sub base_transcript_variation {
    my ($self, $btv) = @_;
    assert_ref($btv, 'Bio::EnsEMBL::Variation::BaseTranscriptVariation') if $btv;
    return $self->variation_feature_overlap($btv);
}

=head2 transcript

  Description: Get the associated Transcript
  Returntype : Bio::EnsEMBL::Transcript
  Exceptions : none
  Status     : Stable

=cut

sub transcript {
    my $self = shift;
    return $self->base_transcript_variation->transcript;
}

=head2 base_variation_feature

  Description: Get the associated BaseVariationFeature
  Returntype : Bio::EnsEMBL::Variation::BaseVariationFeature
  Exceptions : none
  Status     : Stable

=cut

sub base_variation_feature {
    my $self = shift;
    return $self->base_transcript_variation->base_variation_feature;
}

1;