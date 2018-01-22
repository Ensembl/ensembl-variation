=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

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

