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


use Bio::EnsEMBL::DBSQL::DBAdaptor;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);

our @ISA = ('Bio::EnsEMBL::DBSQL::DBAdaptor');


sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);

  my ($dnadb) = rearrange(['DNADB'], @_);

  $self->dnadb($dnadb) if($dnadb);


  return $self;
}


sub dnadb {
  my $self = shift;

  if(@_) {
    my $dnadb = shift;
    if(defined($dnadb) && !$dnadb->isa('Bio::EnsEMBL::DBSQL::DBAdaptor')) {
      throw('Bio::EnsEMBL::DBSQL::DBAdaptor argument expected');
    }
    $self->{'dnadb'} = $dnadb;
  }
  return $self->{'dnadb'};
}

sub get_available_adaptors{
    my %pairs = (
		 'Population' => 'Bio::EnsEMBL::Variation::DBSQL::PopulationAdaptor',
		 'Individual' => 'Bio::EnsEMBL::Variation::DBSQL::IndividualAdaptor',
		 'Variation' => 'Bio::EnsEMBL::Variation::DBSQL::VariationAdaptor',
		 'VariationFeature' => 'Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor',
		 'LDFeatureContainer' => 'Bio::EnsEMBL::Variation::DBSQL::LDFeatureContainerAdaptor',
		 'IndividualGenotype' => 'Bio::EnsEMBL::Variation::DBSQL::IndividualGenotypeAdaptor',
		 'PopulationGenotype' => 'Bio::EnsEMBL::Variation::DBSQL::PopulationGenotypeAdaptor',
		 'TranscriptVariation' => 'Bio::EnsEMBL::Variation::DBSQL::TranscriptVariationAdaptor',
		 'VariationGroup' => 'Bio::EnsEMBL::Variation::DBSQL::VariationGroupAdaptor',
		 'AlleleGroup' => 'Bio::EnsEMBL::Variation::DBSQL::AlleleGroupAdaptor',
		 'VariationGroupFeature' => 'Bio::EnsEMBL::Variation::DBSQL::VariationGroupFeatureAdaptor',
		 'MetaCoordContainer' => 'Bio::EnsEMBL::DBSQL::MetaCoordContainer',
		 'Slice' => 'Bio::EnsEMBL::DBSQL::SliceAdaptor',
		 'CoordSystem'   => 'Bio::EnsEMBL::DBSQL::CoordSystemAdaptor',
		 'AssemblyMapper'       => 'Bio::EnsEMBL::DBSQL::AssemblyMapperAdaptor',
		 'Sequence'             => 'Bio::EnsEMBL::DBSQL::SequenceAdaptor',
		 'Transcript'           => 'Bio::EnsEMBL::DBSQL::TranscriptAdaptor'
		 );
    return (\%pairs);
}

1;
