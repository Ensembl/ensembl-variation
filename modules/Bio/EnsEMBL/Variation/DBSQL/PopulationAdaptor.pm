#
# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::PopulationAdaptor
#
# Copyright (c) 2004 Ensembl
#
# You may distribute this module under the same terms as perl itself
#
#

=head1 NAME

Bio::EnsEMBL::Variation::DBSQL::PopulationAdaptor

=head1 SYNOPSIS

  $db = Bio::EnsEMBL::Variation::DBSQL::DBAdaptor->new(...);

  $pa = $db->get_PopulationAdaptor();

  # Get a Population by its internal identifier
  $pop = $pa->fetch_by_dbID(145);

  # fetch a population by its name
  $pop = $pa->fetch_by_name('PACIFIC'); 

  # fetch all sub populations of a population
  foreach $sub_pop (@{$pa->fetch_all_by_super_Population($pop)}) {
    print $sub_pop->name(), " is a sub population of ", $pop->name(), "\n";
  }

  # fetch all super populations
  foreach $super_pop (@{$pa->fetch_all_by_sub_Population($pop)}) {
    print $pop->name(), " is a sub population of ", $super_pop->name(), "\n";
  }


=head1 DESCRIPTION

This adaptor provides database connectivity for Population objects.
Populations may be retrieved from the Ensembl variation database by
several means using this module.

=head1 AUTHOR - Graham McVicker

=head1 CONTACT

Post questions to the Ensembl development list ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::PopulationAdaptor;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use Bio::EnsEMBL::Variation::Population;

our @ISA = ('Bio::EnsEMBL::DBSQL::BaseAdaptor');



=head2 fetch_by_dbID

  Arg [1]    : int $dbID
  Example    : $pop = $pop_adaptor->fetch_by_dbID(5543);
  Description: Retrieves a Population object via its internal identifier.
               If no such population exists undef is returned.
  Returntype : Bio::EnsEMBL::Variation::Population
  Exceptions : throw if dbID arg is not defined
  Caller     : general, IndividualAdaptor

=cut

sub fetch_by_dbID {
  my $self = shift;
  my $dbID = shift;

  throw('dbID argument expected') if(!defined($dbID));

  my $sth = $self->prepare(q{SELECT population_id, name, size, description
                             FROM   population
                             WHERE  population_id = ?});
  $sth->execute($dbID);
  my $result = $self->_objs_from_sth($sth);
  $sth->finish();

  return undef if(!@$result);

  return $result->[0];
}



=head2 fetch_by_name

  Arg [1]    : string $name
  Example    : $pop = $pop_adaptor->fetch_by_name('NUSPAE:Singapore_HDL');
  Description: Retrieves a population object via its name
  Returntype : Bio::EnsEMBL::Variation::Population
  Exceptions : throw if name argument is not defined
  Caller     : general

=cut

sub fetch_by_name {
  my $self = shift;
  my $name = shift;

  throw('name argument expected') if(!defined($name));

  my $sth = $self->prepare(q{SELECT population_id, name, size, description
                             FROM   population
                             WHERE  name = ?});

  $sth->execute($name);

  my $result = $self->_objs_from_sth($sth);

  $sth->finish();

  return undef if(!@$result);

  return $result->[0];
}


=head2 fetch_all_by_super_Population

  Arg [1]    : Bio::EnsEMBL::Variation::Population $pop
  Example    : foreach $sub_pop (@{$pa->fetch_all_by_super_Population($pop)}) {
                 print $sub_pop->name(), "\n";
               }
  Description: Retrieves all sub populations of a provided population.
  Returntype : Bio::EnsEMBL::Population
  Exceptions : throw on bad argument
  Caller     : general

=cut

sub fetch_all_by_super_Population {
  my $self = shift;
  my $pop  = shift;

  if(!ref($pop) || !$pop->isa('Bio::EnsEMBL::Variation::Population')) {
    throw('Bio::EnsEMBL::Variation::Population argument expected');
  }

  if(!$pop->dbID()) {
    warning("Cannot retrieve sub populations for population without dbID");
    return [];
  }

  my $sth = $self->prepare(q{SELECT p.population_id, p.name, p.size,
                                    p.description
                             FROM   population p, population_structure ps
                             WHERE  p.population_id = ps.sub_population_id
                             AND    ps.super_population_id = ?});

  $sth->execute($pop->dbID());

  my $result = $self->_objs_from_sth($sth);

  $sth->finish();

  return $result;
}



=head2 fetch_all_by_sub_Population

  Arg [1]    : Bio::EnsEMBL::Variation::Population $pop
  Example    : foreach $super_pop (@{$pa->fetch_all_by_sub_Population($pop)}) {
                 print $super_pop->name(), "\n";
               }
  Description: Retrieves all super populations for a provided population
  Returntype : Bio::EnsEMBL::Population
  Exceptions : throw on bad argument
  Caller     : general

=cut

sub fetch_all_by_sub_Population {
  my $self = shift;
  my $pop  = shift;

  if(!ref($pop) || !$pop->isa('Bio::EnsEMBL::Variation::Population')) {
    throw('Bio::EnsEMBL::Variation::Population argument expected');
  }

  if(!$pop->dbID()) {
    warning("Cannot retrieve super populations for population without dbID");
    return [];
  }

  my $sth = $self->prepare(q{SELECT p.population_id, p.name, p.size,
                                    p.description
                             FROM   population p, population_structure ps
                             WHERE  p.population_id = ps.super_population_id
                             AND    ps.sub_population_id = ?});

  $sth->execute($pop->dbID());

  my $result = $self->_objs_from_sth($sth);

  $sth->finish();

  return $result;
}


#
# private method, creates population objects from an executed statement handle
# ordering of columns must be consistant
#
sub _objs_from_sth {
  my $self = shift;
  my $sth  = shift;

  my @pops;

  my ($pop_id, $name, $size, $desc);

  $sth->bind_columns(\$pop_id, \$name, \$size, \$desc);

  while($sth->fetch()) {
    push @pops, Bio::EnsEMBL::Variation::Population->new
      (-dbID => $pop_id,
       -ADAPTOR => $self,
       -NAME => $name,
       -DESCRIPTION => $desc,
       -SIZE => $size);
  }

  return \@pops;
}


1;
