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

# Ensembl module for Bio::EnsEMBL::Variation::Sample
#
# Copyright (c) 2005 Ensembl
#


=head1 NAME

Bio::EnsEMBL::Variation::Sample - An abstract base class to represent
Population, Individual or Strain


=head1 SYNOPSIS

Abstract class - should not be instantiated.  Implementation of
abstract methods must be performed by subclasses.

=head1 DESCRIPTION

This is a base class representing population, individual and strain. This base
class is simply a way of merging similar concepts that should have the same ID


=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::Sample;

use Bio::EnsEMBL::Storable;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw);

our @ISA = ('Bio::EnsEMBL::Storable');


=head2 name

  Arg [1]    : string $newval (optional)
               The new value to set the name attribute to
  Example    : $name = $obj->name()
  Description: Getter/Setter for the name attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub name{
  my $self = shift;
  return $self->{'name'} = shift if(@_);
  return $self->{'name'};
}



=head2 description

  Arg [1]    : string $newval (optional) 
               The new value to set the description attribute to
  Example    : $description = $obj->description()
  Description: Getter/Setter for the description attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub description{
  my $self = shift;
  return $self->{'description'} = shift if(@_);
  return $self->{'description'};
}



=head2 size

  Arg [1]    : int $newval (optional) 
               The new value to set the size attribute to
  Example    : $size = $obj->size()
  Description: Getter/Setter for the size attribute
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub size{
  my $self = shift;
  return $self->{'size'} = shift if(@_);
  return $self->{'size'};
}


1;
