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

package Bio::EnsEMBL::Variation::RegulatoryFeatureVariation;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::RegulatoryFeatureVariationAllele;
use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(overlap);

use base qw(Bio::EnsEMBL::Variation::RegulationVariation);

sub new {
    my $class = shift;

    my %args = @_;

    # swap a '-regulatory_feature' argument for a '-feature' one for the superclass

    for my $arg (keys %args) {
        if (lc($arg) eq '-regulatory_feature') {
            $args{'-feature'} = delete $args{$arg};
        }
    }
   
    # call the superclass constructor
    my $self = $class->SUPER::new(%args) || return undef;
    
    # rebless the alleles from vfoas to rfvas
    map { bless $_, 'Bio::EnsEMBL::Variation::RegulatoryFeatureVariationAllele' } 
        @{ $self->get_all_RegulatoryFeatureVariationAlleles };
    
    return $self;
}

sub regulatory_feature_stable_id {
    my $self = shift;
    return $self->SUPER::_feature_stable_id(@_);
}

sub regulatory_feature {
    my ($self, $rf) = @_;
    return $self->SUPER::feature($rf, 'RegulatoryFeature');
}

sub add_RegulatoryFeatureVariationAllele {
    my $self = shift;
    return $self->SUPER::add_VariationFeatureOverlapAllele(@_);
}

sub get_reference_RegulatoryFeatureVariationAllele {
    my $self = shift;
    return $self->SUPER::get_reference_VariationFeatureOverlapAllele(@_);
}

sub get_all_alternate_RegulatoryFeatureVariationAlleles {
    my $self = shift;
    return $self->SUPER::get_all_alternate_VariationFeatureOverlapAlleles(@_);
}

sub get_all_RegulatoryFeatureVariationAlleles {
    my $self = shift;
    return $self->SUPER::get_all_VariationFeatureOverlapAlleles(@_);
}

1;
