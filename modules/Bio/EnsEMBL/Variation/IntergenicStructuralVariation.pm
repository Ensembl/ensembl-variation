=head1 LICENSE

 Copyright (c) 1999-2012 The European Bioinformatics Institute and
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

package Bio::EnsEMBL::Variation::IntergenicStructuralVariation;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::IntergenicStructuralVariationAllele;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use base qw(Bio::EnsEMBL::Variation::StructuralVariationOverlap);

sub new {
    my $class = shift;

    my %args = @_;

    for my $arg (keys %args) {
        if (lc($arg) eq '-feature') {
            throw("Intergenic variations do not have an associated feature!");
        }
    }   

    # call the superclass constructor
    my $self = $class->SUPER::new(%args) || return undef;
    
    # rebless the alleles from vfoas to ivas
    map { bless $_, 'Bio::EnsEMBL::Variation::IntergenicStructuralVariationAllele' } 
        @{ $self->get_all_IntergenicStructuralVariationAlleles };
    
    return $self;
}

sub feature {
    my $self = shift;
    warning("Intergenic variants do not have an associated feature!") if @_;
    return undef;
}

sub add_IntergenicStructuralVariationAllele {
    my $self = shift;
    return $self->SUPER::add_StructuralVariationOverlapAllele(@_);
}

sub get_reference_IntergenicStructuralVariationAllele {
    my $self = shift;
    return $self->SUPER::get_reference_StructuralVariationOverlapAllele(@_);
}

sub get_all_alternate_IntergenicStructuralVariationAlleles {
    my $self = shift;
    return $self->SUPER::get_all_alternate_StructuralVariationOverlapAlleles(@_);
}

sub get_all_IntergenicStructuralVariationAlleles {
    my $self = shift;
    return $self->SUPER::get_all_StructuralVariationOverlapAlleles(@_);
}

1;
