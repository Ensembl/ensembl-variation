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

package Bio::EnsEMBL::Variation::StructuralVariationOverlap;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::StructuralVariationOverlapAllele;

use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);

use base qw(Bio::EnsEMBL::Variation::BaseVariationFeatureOverlap);

sub new {

    my $class = shift;
    
    my %args = @_;

    # swap a '-structural_variation_feature' argument for a '-base_variation_feature' one for the superclass

    for my $arg (keys %args) {
        if (lc($arg) eq '-structural_variation_feature') {
            $args{'-base_variation_feature'} = delete $args{$arg};
        }
    }

    # call the superclass constructor
    my $self = $class->SUPER::new(%args);
    
    # construct a fake 'allele'
    
    $self->add_StructuralVariationOverlapAllele(
        Bio::EnsEMBL::Variation::StructuralVariationOverlapAllele->new_fast({
            structural_variation_overlap => $self,
            allele_number                => 1,  
        })
    );

    return $self;
}

sub new_fast {
    my ($class, $hashref) = @_;
    
    # swap a 'structural_variation_feature' argument for a 'base_variation_feature' one for the superclass
    
    if ($hashref->{structural_variation_feature}) {
        $hashref->{base_variation_feature} = delete $hashref->{structural_variation_feature};
    }
    
    # and call the superclass

    my $self = $class->SUPER::new_fast($hashref);

#    for my $ssv (@{ $self->structural_variation_feature->structural_variation->get_all_SupportingStructuralVariants }) {
#        for my $ssvf (@{ $ssv->get_all_StructuralVariationFeatures }) {
#            $self->add_StructuralVariationOverlapAllele(
#                Bio::EnsEMBL::Variation::StructuralVariationOverlapAllele->new_fast({
#                    structural_variation_overlap    => $self,
#                })
#            );
#        }
#    }

    unless (@{ $self->get_all_alternate_StructuralVariationOverlapAlleles }) {

        # construct a fake 'allele'
        
        $self->add_StructuralVariationOverlapAllele(
            Bio::EnsEMBL::Variation::StructuralVariationOverlapAllele->new_fast({
                structural_variation_overlap    => $self,
            })
        );
    }

    return $self;
}

sub structural_variation_feature {
    my $self = shift;
    return $self->base_variation_feature(@_);
}

=head2 add_StructuralVariationOverlapAllele

  Arg [1]    : A Bio::EnsEMBL::Variation::StructuralVariationOverlapAllele instance
  Description: Add an allele to this StructuralVariationOverlap
  Returntype : none
  Exceptions : throws if the argument is not the expected type
  Status     : Stable

=cut

sub add_StructuralVariationOverlapAllele {
    my ($self, $svoa) = @_;
    assert_ref($svoa, 'Bio::EnsEMBL::Variation::StructuralVariationOverlapAllele');
    return $self->SUPER::add_BaseVariationFeatureOverlapAllele($svoa);
}

=head2 get_reference_StructuralVariationOverlapAllele

  Description: Get the object representing the reference allele of this StructuralVariationOverlap
  Returntype : Bio::EnsEMBL::Variation::StructuralVariationOverlapAllele instance
  Exceptions : none
  Status     : Stable

=cut

sub get_reference_StructuralVariationOverlapAllele {
    my $self = shift;
    return $self->SUPER::get_reference_BaseVariationFeatureOverlapAllele(@_);
}

=head2 get_all_alternate_StructuralVariationOverlapAlleles

  Description: Get a list of the alternate alleles of this StructuralVariationOverlap
  Returntype : listref of Bio::EnsEMBL::Variation::StructuralVariationOverlapAllele objects
  Exceptions : none
  Status     : Stable

=cut

sub get_all_alternate_StructuralVariationOverlapAlleles {
    my $self = shift;
    return $self->SUPER::get_all_alternate_BaseVariationFeatureOverlapAlleles(@_);
}

=head2 get_all_StructuralVariationOverlapAlleles

  Description: Get a list of the all the alleles, both reference and alternate, of 
               this StructuralVariationOverlap
  Returntype : listref of Bio::EnsEMBL::Variation::StructuralVariationOverlapAllele objects
  Exceptions : none
  Status     : Stable

=cut

sub get_all_StructuralVariationOverlapAlleles {
    my $self = shift;
    return $self->SUPER::get_all_BaseVariationFeatureOverlapAlleles(@_);
}

1;

