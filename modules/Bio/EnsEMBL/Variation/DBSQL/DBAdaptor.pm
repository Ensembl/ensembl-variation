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


our @ISA = ('Bio::EnsEMBL::DBSQL::DBConnection');


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



1;
