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

package Bio::EnsEMBL::Variation::IntergenicVariation;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::IntergenicVariationAllele;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use base qw(Bio::EnsEMBL::Variation::VariationFeatureOverlap);

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
    map { bless $_, 'Bio::EnsEMBL::Variation::IntergenicVariationAllele' } 
        @{ $self->get_all_IntergenicVariationAlleles };
    
    return $self;
}

sub feature {
    my $self = shift;
    warning("Intergenic variants do not have an associated feature!") if @_;
    return undef;
}

sub add_IntergenicVariationAllele {
    my $self = shift;
    return $self->SUPER::add_VariationFeatureOverlapAllele(@_);
}

sub get_reference_IntergenicVariationAllele {
    my $self = shift;
    return $self->SUPER::get_reference_VariationFeatureOverlapAllele(@_);
}

sub get_all_alternate_IntergenicVariationAlleles {
    my $self = shift;
    return $self->SUPER::get_all_alternate_VariationFeatureOverlapAlleles(@_);
}

sub get_all_IntergenicVariationAlleles {
    my $self = shift;
    return $self->SUPER::get_all_VariationFeatureOverlapAlleles(@_);
}

1;
