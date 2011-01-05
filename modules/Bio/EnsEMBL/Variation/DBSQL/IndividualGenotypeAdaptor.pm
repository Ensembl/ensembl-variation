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

#
# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::IndividualGenotypeAdaptor
#
# Copyright (c) 2004 Ensembl
#
# You may distribute this module under the same terms as perl itself
#
#

=head1 NAME

Bio::EnsEMBL::Variation::DBSQL::IndividualGenotypeAdaptor

=head1 SYNOPSIS
  $reg = 'Bio::EnsEMBL::Registry';
  
  $reg->load_registry_from_db(-host => 'ensembldb.ensembl.org',-user => 'anonymous');
  
  $iga = $reg->get_adaptor("human","variation","individualgenotype");
  $sa = $reg->get_adaptor("human","core","slice");
  
  $slice = $sa->fetch_by_region("chromosome", 6, 133088927, 133089926);
  
  foreach my $ig(@{$iga->fetch_all_by_Slice($slice}) {
    print $ig->allele1, "|", $ig->allele2;
  }
  


=head1 DESCRIPTION

This adaptor provides database connectivity for IndividualGenotype objects.
IndividualGenotypes may be retrieved from the Ensembl variation database by
several means using this module.

=head1 METHODS

=cut

package Bio::EnsEMBL::Variation::DBSQL::IndividualGenotypeAdaptor;

use strict;
use warnings;

use vars qw(@ISA);

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Variation::DBSQL::BaseGenotypeAdaptor;

@ISA = ('Bio::EnsEMBL::Variation::DBSQL::BaseGenotypeAdaptor');


=head2 fetch_all_by_Variation

  Arg [1]    : Bio::EnsEMBL::Variation $variation
  Example    : my $var = $variation_adaptor->fetch_by_name( "rs1121" )
               $igtypes = $igtype_adaptor->fetch_all_by_Variation( $var )
  Description: Retrieves a list of individual genotypes for the given Variation.
               If none are available an empty listref is returned.
  Returntype : listref Bio::EnsEMBL::Variation::IndividualGenotype 
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut


sub fetch_all_by_Variation {
    my $self = shift;
    my $variation = shift;
    $self->_multiple(0); #to return data from single and multiple table
    return $self->SUPER::fetch_all_by_Variation($variation);

}

=head2 fetch_all_by_Slice

  Arg [1]    : Bio::EnsEMBL::Slice $slice
  Example    : $igtypes = $igtype_adaptor->fetch_all_by_Slice( $slice )
  Description: Retrieves a list of individual genotypes in a given Slice.
               If none are available an empty listref is returned.
  Returntype : listref Bio::EnsEMBL::Variation::IndividualGenotype 
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub fetch_all_by_Slice{
    my $self = shift;
    my $slice = shift;
    
    $self->_multiple(0);
    my @final_genotypes;
    push @final_genotypes, @{$self->SUPER::fetch_all_by_Slice($slice)};
    $self->_multiple(1);
    $self->SUPER::fetch_all_by_Slice($slice);
    push @final_genotypes, @{$self->SUPER::fetch_all_by_Slice($slice)};
    return \@final_genotypes;
}

1;
