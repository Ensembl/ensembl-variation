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

    # Structural variation feature representing a CNV
    $ssv = Bio::EnsEMBL::Variation::StructuralVariation->new
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

use Bio::EnsEMBL::Storable;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument  qw(rearrange);
use Bio::EnsEMBL::Variation::StructuralVariation;

our @ISA = ('Bio::EnsEMBL::Storable');


=head2 new

  Arg [-dbID] :
    see superclass constructor

  Arg [-ADAPTOR] :
    see superclass constructor

  Arg [-NAME] :
    see superclass constructor
    
  Arg [-STRUCTURAL_VARIATION_ID] :
    see superclass constructor

  Example    :
		
    $svv = Bio::EnsEMBL::Variation::SupportingStructuralVariation->new
       (-name => 'esv25480',
        -structural_variation_id   => $structural_variation->dbID);

  Description: Constructor. Instantiates a new StructuralVariation object.
  Returntype : Bio::EnsEMBL::Variation::StructuralVariation
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);
  my ($dbID, $name, $structural_variation_id) = rearrange([qw(dbID NAME STRUCTURAL_VARIATION_ID)], @_);

  $self->{'dbID'} = $dbID;
  $self->{'name'} = $name;
  $self->{'structural_variation_id'} = $structural_variation_id;
  
  return $self;
}


=head2 name

  Arg [1]    : string $newval (optional)
               The new value to set the name attribute to
  Example    : $name = $obj->name()
  Description: Getter/Setter for the name attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub name{
  my $self = shift;
  return $self->{'name'} = shift if(@_);
  return $self->{'name'};
}


=head2 get_StructuralVariation
  Example    : $sv = $obj->get_StructuralVariation()
  Description: Getter of the structural variation supported by the supporting evidence. 
  Returntype : Bio::EnsEMBL::Variation::StructuralVariation
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub get_StructuralVariation{
  my $self = shift;

	my $sva = $self->{'adaptor'}->db()->get_StructuralVariationAdaptor();
	return $sva->fetch_by_dbID($self->{'structural_variation_id'});
}

=head2 is_structural_variation
  Example    : $sv = $obj->is_structural_variation()
  Description: Getter to determine if the supporting evidence is also a structural variant 
  Returntype : Bio::EnsEMBL::Variation::StructuralVariation
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub is_structural_variation{
  my $self = shift;

	my $sva = $self->{'adaptor'}->db()->get_StructuralVariationAdaptor();
	return $sva->fetch_by_name($self->{'name'});
}
1;

