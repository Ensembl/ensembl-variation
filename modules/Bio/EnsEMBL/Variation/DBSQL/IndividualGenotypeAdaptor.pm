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
  $ia = $reg->get_adaptor("human","variation","individual");

  # Get an IndividualGenotype by its internal identifier
  $igtype = $ia->fetch_by_dbID(145);

  # Get all individual genotypes for an individual
  $ind = $ia->fetch_all_by_Individual(1219);

  foreach $igtype (@{$iga->fetch_all_by_Individual($ind)}) {
    print $igtype->variation()->name(),  ' ',
          $igtype->allele1(), '/', $igtype->allele2(), "\n";
  }



=head1 DESCRIPTION

This adaptor provides database connectivity for IndividualGenotype objects.
IndividualGenotypes may be retrieved from the Ensembl variation database by
several means using this module.

=head1 AUTHOR - Graham McVicker

=head1 CONTACT

Post questions to the Ensembl development list ensembl-dev@ebi.ac.uk

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

=head2 fetch_all_by_Individual

  Arg [1]    : Bio::EnsEMBL::Variation::Individual
  Example    : $ind = $ind_adaptor->fetch_by_dbID(1345);
               @gtys = $igty_adaptor->fetch_all_by_Individual($ind);
  Description: Retrieves all genotypes which are stored for a specified
               individual.
  Returntype : Bio::EnsEMBL::Variation::IndividualGenotype
  Exceptions : throw on incorrect argument
  Caller     : general
  Status     : At Risk

=cut

sub fetch_all_by_Individual {
  my $self = shift;
  my $ind = shift;
  

  $self->_multiple(0); #to return data from single and multiple genotype table
  return $self->SUPER::fetch_all_by_Individual($ind);
}


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
