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

package Bio::EnsEMBL::Variation::ExternalFeatureVariation;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::ExternalFeatureVariationAllele;
use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(overlap);

use base qw(Bio::EnsEMBL::Variation::RegulationVariation);

sub new {
    my $class = shift;

    my %args = @_;

    # swap a '-external_feature' argument for a '-feature' one for the superclass

    for my $arg (keys %args) {
        if (lc($arg) eq '-external_feature') {
            $args{'-feature'} = delete $args{$arg};
        }
    }   

    # call the superclass constructor
    my $self = $class->SUPER::new(%args) || return undef;
    
    # rebless the alleles from vfoas to efvas
    map { bless $_, 'Bio::EnsEMBL::Variation::ExternalFeatureVariationAllele' } 
        @{ $self->get_all_ExternalFeatureVariationAlleles };
    
    return $self;
}

sub external_feature_stable_id {
    my $self = shift;
    return $self->SUPER::feature_stable_id(@_);
}

sub external_feature {
    my ($self, $ef) = @_;
    return $self->SUPER::feature($ef, 'ExternalFeature');
}

sub add_ExternalFeatureVariationAllele {
    my $self = shift;
    return $self->SUPER::add_VariationFeatureOverlapAllele(@_);
}

sub get_reference_ExternalFeatureVariationAllele {
    my $self = shift;
    return $self->SUPER::get_reference_VariationFeatureOverlapAllele(@_);
}

sub get_all_alternate_ExternalFeatureVariationAlleles {
    my $self = shift;
    return $self->SUPER::get_all_alternate_VariationFeatureOverlapAlleles(@_);
}

sub get_all_ExternalFeatureVariationAlleles {
    my $self = shift;
    return $self->SUPER::get_all_VariationFeatureOverlapAlleles(@_);
}

sub target_feature_stable_id {
    my ($self, $target_feature_stable_id) = @_;
    
    $self->{target_feature_stable_id} = $target_feature_stable_id if $target_feature_stable_id;
    
    unless ($self->{target_feature_stable_id}) {
        
        # try to fetch it from the funcgen database

        if (my $species = $self->{adaptor}->db->species) {
            for my $entry (
                @{ $self->external_feature->get_all_DBEntries($species.'_core_Transcript') }, 
                @{ $self->external_feature->get_all_DBEntries($species.'_core_Gene') } ) {
                if (my $id = $entry->primary_id) {

                    $self->{target_feature_stable_id} = $id;

                    # there should never be more than one, so we last out of the loop
                    
                    last;
                }
            }
        }
        else {
            warn "Failed to get species from adaptor";
        }
    }
    
    return $self->{target_feature_stable_id};
}

1;
