#
# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::DBAdaptor
#
# Copyright (c) 2004 Ensembl
#
# You may distribute this module under the same terms as perl itself
#
#

=head1 NAME

Bio::EnsEMBL::DBSQL::DBAdaptor

=head1 SYNOPSIS



=head1 DESCRIPTION

This module provides a connection to an Ensembl variation database and
provides a means to obtain ObjectAdaptors.

=head1 AUTHOR - Graham McVicker

=head1 CONTACT

Post questions to the Ensembl development list ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;


package Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;


use Bio::EnsEMBL::DBSQL::DBConnection;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);

our @ISA = ('Bio::EnsEMBL::DBSQL::DBConnection');


sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);

  my ($dnadb) = rearrange(['DNADB'], @_);

  $self->dnadb($dnadb) if($dnadb);


  return $self;
}


sub dnadb {
  my $self = shift;

  $self = $self->_obj() if($self->isa('Bio::EnsEMBL::Container'));

  if(@_) {
    my $dnadb = shift;
    if(defined($dnadb) && !$dnadb->isa('Bio::EnsEMBL::DBSQL::DBAdaptor')) {
      throw('Bio::EnsEMBL::DBSQL::DBAdaptor argument expected');
    }
    $self->{'dnadb'} = $dnadb;
  }
  return $self->{'dnadb'};
}


sub get_PopulationAdaptor {
  my $self = shift;
  return $self->_get_adaptor
    ('Bio::EnsEMBL::Variation::DBSQL::PopulationAdaptor');
}

sub get_IndividualAdaptor {
  my $self = shift;
  return $self->_get_adaptor
    ('Bio::EnsEMBL::Variation::DBSQL::IndividualAdaptor');
}


sub get_VariationAdaptor {
  my $self = shift;
  return $self->_get_adaptor
    ('Bio::EnsEMBL::Variation::DBSQL::VariationAdaptor');
}

sub get_VariationFeatureAdaptor {
  my $self = shift;
  return $self->_get_adaptor
    ('Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor');
}

sub get_IndividualGenotypeAdaptor {
  my $self = shift;
  return $self->_get_adaptor
    ('Bio::EnsEMBL::Variation::DBSQL::IndividualGenotypeAdaptor');
}

sub get_PopulationGenotypeAdaptor {
  my $self = shift;
  return $self->_get_adaptor
    ('Bio::EnsEMBL::Variation::DBSQL::PopulationGenotypeAdaptor');
}


sub get_TranscriptVariationAdaptor {
  my $self = shift;
  return $self->_get_adaptor
    ('Bio::EnsEMBL::Variation::DBSQL::TranscriptVariationAdaptor');
}


sub get_SliceAdaptor {
  my $self = shift;
  if(!$self->dnadb()) {
    throw('Cannot obtain SliceAdaptor without attached dnadb');
  }
  return $self->dnadb->get_SliceAdaptor();
}

sub get_CoordSystemAdaptor {
  my $self = shift;
  if(!$self->dnadb()) {
    throw('Cannot obtain CoordSystemAdaptor without attached dnadb');
  }
  return $self->dnadb->get_CoordSystemAdaptor();
}


sub get_AssemblyMapperAdaptor {
  my $self = shift;
  if(!$self->dnadb()) {
    throw('Cannot obtain AssemblyMapperAdaptor without attached dnadb');
  }
  return $self->dnadb->get_AssemblyMapperAdaptor();
}


sub get_SequenceAdaptor {
  my $self = shift;
  if(!$self->dnadb()) {
    throw("Cannot obtain SequenceAdaptor without attached dnadb");
  }
  return $self->dnadb->get_SequenceAdaptor();
}


sub get_TranscriptAdaptor {
  my $self = shift;
  if(!$self->dnadb()) {
    throw("Cannot obtain TranscriptAdaptor without attached dnadb");
  }
  return $self->dnadb->get_TranscriptAdaptor();
}


sub get_MetaCoordContainer {
  my $self = shift;
  return $self->_get_adaptor('Bio::EnsEMBL::DBSQL::MetaCoordContainer');
}

1;
