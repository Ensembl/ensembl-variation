=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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
 developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

 Questions may also be sent to the Ensembl help desk at
 <http://www.ensembl.org/Help/Contact>.

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

use base qw(Bio::EnsEMBL::Variation::BaseVariationFeatureOverlapAllele);

=head2 base_transcript_variation

  Description: Get/set the associated BaseTranscriptVariation
  Returntype : Bio::EnsEMBL::Variation::BaseTranscriptVariation
  Exceptions : throws if the argument is the wrong type
  Status     : Stable

=cut

sub base_transcript_variation {
    my ($self, $btv) = @_;
    assert_ref($btv, 'Bio::EnsEMBL::Variation::BaseTranscriptVariation') if $btv;
    return $self->base_variation_feature_overlap($btv);
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