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

# Ensembl module for Bio::EnsEMBL::Variation::SupportingStructuralVariation
#
# Copyright (c) 2011 Ensembl
#


=head1 NAME

Bio::EnsEMBL::Variation::SupportingStructuralVariation - A supporting evidence for a structural variation.

=head1 SYNOPSIS

    # Supporting evidence of a structural variation
    $ssv = Bio::EnsEMBL::Variation::SupportingStructuralVariation->new
       (-supporting_structural_evidence => 'ssv001'
        -structural_variation   => $structural_variation);

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
use Bio::EnsEMBL::Utils::Exception qw(throw warning deprecate);

our @ISA = ('Bio::EnsEMBL::Variation::BaseStructuralVariation');


sub new {
	my $caller = shift;
	my $class = ref($caller) || $caller;
	
	my $self = Bio::EnsEMBL::Variation::BaseStructuralVariation->new(@_);
	return(bless($self, $class));
}

=head2 name

  Arg [1]    : string $newval (optional)
               The new value to set the name attribute to
  Example    : $name = $obj->name()
  Description: Getter/Setter for the name attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : DEPRECATED: use the variation_name method

=cut

sub name{
  my $self = shift;
	deprecate('Use the method "variation_name" instead');
  return $self->{'variation_name'} = shift if(@_);
  return $self->{'variation_name'};
}


=head2 get_StructuralVariation
  Example    : $ssv = $obj->get_StructuralVariation()
  Description: Getter of the structural variation supported by the supporting evidence. 
  Returntype : Bio::EnsEMBL::Variation::StructuralVariation
  Exceptions : none
  Caller     : general
  Status     : DEPRECATED: use the get_all_StructuralVariations method

=cut

sub get_StructuralVariation {
  my $self = shift;
	deprecate('Use the method "get_all_StructuralVariations" instead');
}


=head2 get_all_StructuralVariations
  Example    : $ssv = $obj->get_all_StructuralVariations()
  Description: Getter of the structural variations supported by the supporting evidence. 
  Returntype : reference to list of Bio::EnsEMBL::Variation::StructuralVariation objects
  Exceptions : none
  Caller     : general
  Status     : At Risk

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

=head2 is_structural_variation
  Example    : $sv = $obj->is_structural_variation()
  Description: Getter to determine if the supporting evidence is also a structural variant 
  Returntype : Bio::EnsEMBL::Variation::StructuralVariation
  Exceptions : none
  Caller     : general
  Status     : DEPRECATED: no more used

=cut

sub is_structural_variation{
  my $self = shift;
	deprecate('Method no more used');
}
1;

