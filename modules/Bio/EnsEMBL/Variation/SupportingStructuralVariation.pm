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

