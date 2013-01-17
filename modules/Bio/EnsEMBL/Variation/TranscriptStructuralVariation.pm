=head1 LICENSE

 Copyright (c) 1999-2013 The European Bioinformatics Institute and
 Genome Research Limited.  All rights reserved.

 This software is distributed under a modified Apache license.
 For license details, please see

   http://www.ensembl.org/info/about/legal/code_licence.html

=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <dev@ensembl.org>.

 Questions may also be sent to the Ensembl help desk at
 <helpdesk@ensembl.org>.

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

