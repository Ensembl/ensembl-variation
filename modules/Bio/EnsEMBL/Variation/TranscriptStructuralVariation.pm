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

Bio::EnsEMBL::Variation::TranscriptStructuralVariation

=head1 SYNOPSIS

    use Bio::EnsEMBL::Variation::TranscriptStructuralVariation;

    my $tv = Bio::EnsEMBL::Variation::TranscriptVariation->new(
        -transcript                   => $transcript,
        -structural_variation_feature => $svf
    );

=head1 DESCRIPTION

A TranscriptStructuralVariation object represents a structural variation feature
which is in close proximity to an Ensembl transcript. A
TranscriptStructuralVariation object has several attributes which define the
relationship of the variation to the transcript.

=cut

package Bio::EnsEMBL::Variation::TranscriptStructuralVariation;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::TranscriptStructuralVariationAllele;

use base qw(Bio::EnsEMBL::Variation::StructuralVariationOverlap Bio::EnsEMBL::Variation::BaseTranscriptVariation);

sub new {

    my $class = shift;
    
    my %args = @_;

    # swap a '-transcript' argument for a '-feature' one for the superclass

    for my $arg (keys %args) {
        if (lc($arg) eq '-transcript') {
            $args{'-feature'} = delete $args{$arg};
        }
    }

    # call the superclass constructor
    my $self = $class->SUPER::new(%args) || return undef;

    # rebless the alleles from vfoas to tvas
    map { bless $_, 'Bio::EnsEMBL::Variation::TranscriptStructuralVariationAllele' } 
        @{ $self->get_all_TranscriptStructuralVariationAlleles };
    
    return $self;
}

=head2 add_TranscriptStructuralVariationAllele

  Arg [1]    : A Bio::EnsEMBL::Variation::TranscriptStructuralVariationAllele instance
  Description: Add an allele to this TranscriptStructuralVariation
  Returntype : none
  Exceptions : throws if the argument is not the expected type
  Status     : Stable

=cut

sub add_TranscriptStructuralVariationAllele {
    my ($self, $svoa) = @_;
    assert_ref($svoa, 'Bio::EnsEMBL::Variation::TranscriptStructuralVariationAllele');
    return $self->SUPER::add_BaseVariationFeatureOverlapAllele($svoa);
}

=head2 get_reference_TranscriptStructuralVariationAllele

  Description: Get the object representing the reference allele of this TranscriptStructuralVariation
  Returntype : Bio::EnsEMBL::Variation::TranscriptStructuralVariationAllele instance
  Exceptions : none
  Status     : Stable

=cut

sub get_reference_TranscriptStructuralVariationAllele {
    my $self = shift;
    return $self->SUPER::get_reference_BaseVariationFeatureOverlapAllele(@_);
}

=head2 get_all_alternate_TranscriptStructuralVariationAlleles

  Description: Get a list of the alternate alleles of this TranscriptStructuralVariation
  Returntype : listref of Bio::EnsEMBL::Variation::TranscriptStructuralVariationAllele objects
  Exceptions : none
  Status     : Stable

=cut

sub get_all_alternate_TranscriptStructuralVariationAlleles {
    my $self = shift;
    return $self->SUPER::get_all_alternate_BaseVariationFeatureOverlapAlleles(@_);
}

=head2 get_all_TranscriptStructuralVariationAlleles

  Description: Get a list of the all the alleles, both reference and alternate, of 
               this TranscriptStructuralVariation
  Returntype : listref of Bio::EnsEMBL::Variation::TranscriptStructuralVariationAllele objects
  Exceptions : none
  Status     : Stable

=cut

sub get_all_TranscriptStructuralVariationAlleles {
    my $self = shift;
    return $self->SUPER::get_all_BaseVariationFeatureOverlapAlleles(@_);
}

1;

