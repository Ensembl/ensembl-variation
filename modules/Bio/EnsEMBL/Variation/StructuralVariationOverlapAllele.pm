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

# required for intron overlap checking
sub _get_differing_regions {
  my $self = shift;
  my $bvf = $self->base_variation_feature;
  return $self->{_differing_regions} ||= [{ s => 0, e => ($bvf->{end} - $bvf->{start}) }];
}

1;

