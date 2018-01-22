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
    unless($args{'-feature'} ||= delete $args{'-regulatory_feature'}) {
        for my $arg (keys %args) {
            if (lc($arg) eq '-regulatory_feature') {
                $args{'-feature'} = delete $args{$arg};
            }
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
