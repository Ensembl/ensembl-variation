#
# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::PopulationGenotypeAdaptor
#
# Copyright (c) 2004 Ensembl
#
# You may distribute this module under the same terms as perl itself
#
#

=head1 NAME

Bio::EnsEMBL::Variation::DBSQL::PopulationGenotypeAdaptor

=head1 SYNOPSIS

  $db = Bio::EnsEMBL::Variation::DBSQL::DBAdaptor->new(...);

  $pga = $db->get_PopulationGenotypeAdaptor();
  $pa  = $db->get_PopulationAdaptor();

  # Get a PopulationGenotype by its internal identifier
  $pgtype = $ia->fetch_by_dbID(145);

  print $pgtype->population->name(), " ",
        $pgtype->allele1(), ' ', $pgtype->allele2(), ' ', $pgtype->frequency();

  # Get all population genotypes for an population
  $pop = $pa->fetch_by_dbID(1219);

  foreach $pgtype (@{$pga->fetch_all_by_Population($pop)}) {
    print $pgtype->variation()->name(),  ' ',
          $pgtype->frequency();
          $pgtype->allele1(), '/', $pgtype->allele2(), "\n";
  }



=head1 DESCRIPTION

This adaptor provides database connectivity for PopulationGenotype objects.
PopulationGenotypes may be retrieved from the Ensembl variation database by
several means using this module.

=head1 AUTHOR - Graham McVicker

=head1 CONTACT

Post questions to the Ensembl development list ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::PopulationGenotypeAdaptor;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use Bio::EnsEMBL::Variation::PopulationGenotype;

our @ISA = ('Bio::EnsEMBL::DBSQL::BaseAdaptor');



=head2 fetch_by_dbID

  Arg [1]    : int $dbID
  Example    : $pgtype = $pgtype_adaptor->fetch_by_dbID(15767);
  Description: Retrieves a population genotype via its unique internal
               identifier.  undef is returned if no such population genotype
               exists.
  Returntype : Bio::EnsEMBL::Variation::Variation::PopulationGenotype or undef
  Exceptions : throw if no dbID argument is provided
  Caller     : general

=cut

sub fetch_by_dbID {
  my $self = shift;
  my $dbID = shift;

  my $sth = $self->prepare
    (q{SELECT variation_id, allele_1, allele_2, frequency, population_id
       FROM   population_genotype
       WHERE  population_genotype_id = ?});

  $sth->execute($dbID);

  my $row = $sth->fetchrow_arrayref();
  $sth->finish();

  return undef if(!$row);

  my ($var_id, $allele1, $allele2, $freq, $pop_id) = @$row;

  my $va = $self->db()->get_VariationAdaptor();
  my $var = $va->fetch_by_dbID($var_id);

  my $pa = $self->db()->get_PopulationAdaptor();
  my $pop = $pa->fetch_by_dbID($pop_id);

  return Bio::EnsEMBL::Variation::PopulationGenotype->new
    (-dbID    => $dbID,
     -adaptor => $self,
     -allele1 => $allele1,
     -allele2 => $allele2,
     -variation => $var,
     -frequency => $freq,
     -population => $pop);
}



=head2 fetch_all_by_Population

  Arg [1]    : Bio::EnsEMBL::Variation::Population
  Example    : $pop = $pop_adaptor->fetch_by_dbID(1345);
               @gtys = $pgty_adaptor->fetch_all_by_Population($pop);
  Description: Retrieves all genotypes which are stored for a specified
               population.
  Returntype : Bio::EnsEMBL::Variation::PopulationGenotype
  Exceptions : throw on incorrect argument
  Caller     : general

=cut

sub fetch_all_by_Population {
  my $self = shift;
  my $pop = shift;

  if(!ref($pop) || !$pop->isa('Bio::EnsEMBL::Variation::Population')) {
    throw('Bio::EnsEMBL::Variation::Population argument expected');
  }

  if(!defined($pop->dbID())) {
    warning("Cannot retrieve genotypes for population without set dbID");
    return [];
  }

  my $sth = $self->prepare
    (q{SELECT population_genotype_id, variation_id, allele_1, allele_2,
              frequency
       FROM   population_genotype
       WHERE  population_id = ?});

  $sth->execute($pop->dbID());

  my %variation_hash;

  my ($pgtype_id, $var_id, $allele1, $allele2, $freq);
  $sth->bind_columns(\$pgtype_id, \$var_id, \$allele1, \$allele2, \$freq);

  my @results;

  while($sth->fetch()) {
    my $pgtype = Bio::EnsEMBL::Variation::PopulationGenotype->new
      (-dbID => $pgtype_id,
       -adaptor => $self,
       -allele1 => $allele1,
       -allele2 => $allele2,
       -frequency => $freq,
       -population => $pop);
    $variation_hash{$var_id} ||= [];
    push @{$variation_hash{$var_id}}, $pgtype;
    push @results, $pgtype;
  }

  # get all variations in one query (faster)
  # and add to already created genotypes
  my @var_ids = keys %variation_hash;
  my $va = $self->db()->get_VariationAdaptor();
  my $vars = $va->fetch_all_by_dbID_list(\@var_ids);

  foreach my $v (@$vars) {
    foreach my $pgty (@{$variation_hash{$v->dbID()}}) {
      $pgty->variation($v);
    }
  }

  return \@results;
}

1;
