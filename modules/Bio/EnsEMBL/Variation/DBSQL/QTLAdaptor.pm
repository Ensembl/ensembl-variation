=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut


=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

 Questions may also be sent to the Ensembl help desk at
 <http://www.ensembl.org/Help/Contact>.

=cut


# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::PhenotypeFeatureAdaptor
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