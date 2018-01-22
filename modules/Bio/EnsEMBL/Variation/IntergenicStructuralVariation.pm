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
