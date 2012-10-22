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

# Ensembl module for Bio::EnsEMBL::Variation::Phenotype
#
# Copyright (c) 2011 Ensembl
#


=head1 NAME

Bio::EnsEMBL::Variation::Phenotype - Ensembl representation of a phenotype.

=head1 SYNOPSIS

    $study = Bio::EnsEMBL::Variation::Study->new(-DESCRIPTION => 'Hemostatic factors and hematological phenotypes');

=head1 DESCRIPTION

This is a class representing a phenotype from the ensembl-variation database.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::Phenotype;

use Bio::EnsEMBL::Storable;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument  qw(rearrange);

our @ISA = ('Bio::EnsEMBL::Storable');


=head2 new
    Arg [-DESCRIPTION] :
    phenotype description
			
  Example    :
		
    $phenotype = Bio::EnsEMBL::Variation::Phenotype->new(-DESCRIPTION => 'Hemostatic factors and hematological phenotypes');

  Description: Constructor. Instantiates a new Phenotype object.
  Returntype : Bio::EnsEMBL::Variation::Phenotype
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub new {
    my $caller = shift;
    my $class  = ref($caller) || $caller;
    my $self = $class->SUPER::new(@_);
    my ($dbID, $description) = rearrange([qw(dbID DESCRIPTION)], @_);
    $self = {
        'dbID'        => $dbID,
        'description' => $description,
    };
    return bless $self, $class;
}

=head2 id

  Example    : $name = $obj->id()
  Description: Getter/Setter for the id attribute
  Returntype : integer
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub id {
    my $self = shift;
    return $self->{'dbID'} = shift if(@_);
    return $self->{'dbID'};
}

=head2 description

  Example    : $name = $obj->description()
  Description: Getter/Setter for the description attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub description {
    my $self = shift;
    return $self->{'description'} = shift if(@_);
    return $self->{'description'};
}

