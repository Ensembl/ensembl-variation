=head1 LICENSE

 Copyright (c) 1999-2011 The European Bioinformatics Institute and
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

package Bio::EnsEMBL::Variation::ExternalFeatureVariation;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::ExternalFeatureVariationAllele;
use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(overlap);

use base qw(Bio::EnsEMBL::Variation::RegulationVariation);

sub new {
    my $class = shift;
    
    # call the superclass constructor
    my $self = $class->SUPER::new(@_) || return undef;
    
    # rebless the alleles from vfoas to efvas
    map { bless $_, 'Bio::EnsEMBL::Variation::ExternalFeatureVariationAllele' } @{ $self->alleles };
    
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

sub target_feature_stable_id {
    my ($self, $target_feature_stable_id) = @_;
    
    $self->{target_feature_stable_id} = $target_feature_stable_id if $target_feature_stable_id;
    
    unless ($self->{target_feature_stable_id}) {
        
        if (my $species = $self->{adaptor}->db->species) {
            for my $entry (
                @{ $self->feature->get_all_DBEntries($species.'_core_Transcript') }, 
                @{ $self->feature->get_all_DBEntries($species.'_core_Gene') } ) {
                if (my $id = $entry->primary_id) {
                    $self->{target_feature_stable_id} = $id;
                    last;
                }
            }
        }
        else {
            warn "No ad"
        }
    }
    
    return $self->{target_feature_stable_id};
}

1;
