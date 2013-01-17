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

#
# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::PhenotypeAdaptor
#
# Copyright (c) 2004 Ensembl
#
# You may distribute this module under the same terms as perl itself
#
#

=head1 NAME

Bio::EnsEMBL::Variation::DBSQL::PhenotypeAdaptor

=head1 SYNOPSIS
  $reg = 'Bio::EnsEMBL::Registry';
  
  $reg->load_registry_from_db(-host => 'ensembldb.ensembl.org',-user => 'anonymous');
  
  $pa = $reg->get_adaptor("human","variation","phenotype");

  # Get a list of all phenotypes.
  $phenotypes = $pa->fetch_all();

=head1 DESCRIPTION

This adaptor provides database connectivity for Phenotype objects.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::PhenotypeAdaptor;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Variation::Phenotype;

our @ISA = ('Bio::EnsEMBL::DBSQL::BaseAdaptor');

=head2 fetch_all
    
    Example         : my $phenotypes = $phenotype_adaptor->fetch_all();
    Description     : Retrieves an array of all phenotyes.
    Returntype      : listref of Bio::EnsEMBL::Variation::Phenotype. 
    Exceptions      : none
    Caller          : general
    Status          : Stable

=cut

sub fetch_all {
    my $self = shift;
    my @phenotypes;
    
    my $sth = $self->prepare(qq{SELECT phenotype_id, description from phenotype});
    $sth->execute();
    my ($dbID, $phenotype_description);
    $sth->bind_columns(\$dbID, \$phenotype_description);
    while ($sth->fetch) {
        push @phenotypes, Bio::EnsEMBL::Variation::Phenotype->new(
        -dbID        => $dbID,
        -DESCRIPTION => $phenotype_description,);
    }    
    
    return \@phenotypes;
}
