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

# Ensembl module for Bio::EnsEMBL::Variation::Phenotype
#
# Copyright (c) 2011 Ensembl
#


=head1 NAME

Bio::EnsEMBL::Variation::Phenotype - Ensembl representation of a phenotype.

=head1 SYNOPSIS

    my $phenotype = Bio::EnsEMBL::Variation::Phenotype->new(-NAME => 'Type I diabetes');

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
    my ($dbID, $description, $name) = rearrange([qw(dbID DESCRIPTION NAME)], @_);
    $self = {
        'dbID'        => $dbID,
        'description' => $description,
        'name'        => $name,
    };
    return bless $self, $class;
}

sub new_fast {
  my $class = shift;
  my $hashref = shift;
  return bless $hashref, $class;
}

=head2 dbID

  Example    : $name = $obj->dbID()
  Description: Getter/Setter for the dbIDattribute
  Returntype : integer
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub dbID {
    my $self = shift;
    return $self->{'dbID'} = shift if(@_);
    return $self->{'dbID'};
}

=head2 name

  Example    : $name = $obj->name()
  Description: Getter/Setter for the name attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub name {
    my $self = shift;
    return $self->{'name'} = shift if(@_);
    return $self->{'name'};
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

