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

# Ensembl module for Bio::EnsEMBL::Variation::SupportingStructuralVariation
#
#


=head1 NAME

Bio::EnsEMBL::Variation::SupportingStructuralVariation - A supporting evidence for a structural variation.

=head1 SYNOPSIS

    # Supporting evidence of a structural variation
    $ssv = Bio::EnsEMBL::Variation::SupportingStructuralVariation->new
       (-supporting_structural_evidence => 'ssv001'
        -structural_variation   => $structural_variation); # a StructuralVariation object

    ...


=head1 DESCRIPTION

This is a class representing the supporting evidence of a structural variation
from the ensembl-variation database.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::SupportingStructuralVariation;

use Bio::EnsEMBL::Variation::BaseStructuralVariation;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

our @ISA = ('Bio::EnsEMBL::Variation::BaseStructuralVariation');


sub new {
	my $caller = shift;
	my $class = ref($caller) || $caller;
	
	my $self = Bio::EnsEMBL::Variation::BaseStructuralVariation->new(@_);
	return(bless($self, $class));
}


=head2 get_all_StructuralVariations
  Example    : $ssv = $obj->get_all_StructuralVariations()
  Description: Getter of the structural variations supported by the supporting evidence. 
  Returntype : reference to list of Bio::EnsEMBL::Variation::StructuralVariation objects
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_StructuralVariations {
  my $self = shift;

	if(defined $self->{'adaptor'}) {
		my $sva = $self->{'adaptor'}->db()->get_StructuralVariationAdaptor();
		return $sva->fetch_all_by_supporting_evidence($self);
	}
	else {
  	warn("No variation database attached");
  	return [];
  }
}

1;

