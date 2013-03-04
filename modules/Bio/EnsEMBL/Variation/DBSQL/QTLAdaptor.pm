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


# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::PhenotypeFeatureAdaptor
#
# Copyright (c) 2004 Ensembl
#
# You may distribute this module under the same terms as perl itself
#
#

=head1 NAME

Bio::EnsEMBL::Variation::DBSQL::QTLAdaptor

=head1 SYNOPSIS
  my $reg = 'Bio::EnsEMBL::Registry';
  my $qtl_adaptor = $reg->get_adaptor("human","variation","qtl");

  my $qtl = $qtl_adaptor->fetch_by_stable_id('QTL00123');

=head1 DESCRIPTION

This adaptor creates placeholder QTL objects.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::QTLAdaptor;

use Bio::EnsEMBL::Variation::QTL;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;

our @ISA = ('Bio::EnsEMBL::DBSQL::BaseAdaptor');

sub fetch_by_stable_id {
	my $self = shift;
	
	return Bio::EnsEMBL::Variation::QTL->new_fast({
		stable_id => shift;
	})
}

1;