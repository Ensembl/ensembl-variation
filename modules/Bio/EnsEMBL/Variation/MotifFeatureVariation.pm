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

package Bio::EnsEMBL::Variation::MotifFeatureVariation;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::MotifFeatureVariationAllele;
use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(overlap);

use base qw(Bio::EnsEMBL::Variation::RegulationVariation);

sub new {
    my $class = shift;

    my %args = @_;

    # swap a '-motif_feature' argument for a '-feature' one for the superclass

    for my $arg (keys %args) {
        if (lc($arg) eq '-motif_feature') {
            $args{'-feature'} = delete $args{$arg};
        }
    }

    # call the superclass constructor
    my $self = $class->SUPER::new(%args) || return undef;
    
    # rebless the alleles from vfoas to mfvas
    map { bless $_, 'Bio::EnsEMBL::Variation::MotifFeatureVariationAllele' } 
        @{ $self->get_all_MotifFeatureVariationAlleles };
    
    return $self;
}

sub motif_feature_stable_id {
    my $self = shift;
    return $self->SUPER::feature_stable_id(@_);
}

sub motif_feature {
    my ($self, $mf) = @_;
    return $self->SUPER::feature($mf, 'MotifFeature');
}

sub add_MotifFeatureVariationAllele {
    my $self = shift;
    return $self->SUPER::add_VariationFeatureOverlapAllele(@_);
}

sub get_reference_MotifFeatureVariationAllele {
    my $self = shift;
    return $self->SUPER::get_reference_VariationFeatureOverlapAllele(@_);
}

sub get_all_alternate_MotifFeatureVariationAlleles {
    my $self = shift;
    return $self->SUPER::get_all_alternate_VariationFeatureOverlapAlleles(@_);
}

sub get_all_MotifFeatureVariationAlleles {
    my $self = shift;
    return $self->SUPER::get_all_VariationFeatureOverlapAlleles(@_);
}

sub _motif_feature_seq {
    my $self = shift;
    
    my $mf = $self->motif_feature;
    
    my $mf_seq = $mf->{_variation_effect_feature_cache}->{seq} ||= $mf->seq;
    
    return $mf_seq;
}

1;
