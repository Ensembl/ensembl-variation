#
# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::IndividualAdaptor
#
# Copyright (c) 2004 Ensembl
#
# You may distribute this module under the same terms as perl itself
#
#

=head1 NAME

Bio::EnsEMBL::DBSQL::IndividualAdaptor

=head1 SYNOPSIS

  my $db = Bio::EnsEMBL::Variation::DBSQL::DBAdaptor->new(...);

  my $ia = $db->get_IndividualAdaptor();

  # Get an individual by its internal identifier
  my $ind = $ia->fetch_by_dbID(52);

  # Get all individuals with a particular name
  foreach my $ind (@{$ia->fetch_all_by_name('PKH053(M)')}) {
    print "Individual ", $ind->name(), "\n";
  }

  # get all individuals from a population
  my $pop = $pop_adaptor->fetch_by_name('PACIFIC');
  foreach my $ind (@{$ia->fetch_all_by_Population($pop)}) {
    print $ind->name(), "\n";
  }

  # get all children of an individual
  foreach my $child (@{$ia->fetch_all_by_parent($ind)}) {
    print $child->name(), " is a child of ", $ind->name(), "\n";
  }


=head1 DESCRIPTION

This adaptor provides database connectivity for Individual objects.
Individuals may be retrieved from the ensembl variation database by
several means using this module.

=head1 AUTHOR - Graham McVicker

=head1 CONTACT

Post questions to the Ensembl development list ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::IndividualAdaptor;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use Bio::EnsEMBL::Variation::Individual;

our @ISA = ('Bio::EnsEMBL::DBSQL::BaseAdaptor');



=head2 fetch_by_dbID

  Arg [1]    : int $dbID - internal identifier of the individual to fetch
  Example    : my $ind = $ind_adptor->fetch_by_dbID(451);
  Description: Retrieves an individual via its internal identifier.
               If no such individual exists in the database undef is returned.
  Returntype : Bio::EnsEMBL::Individual or undef
  Exceptions : throw if no argument provided
  Caller     : general, 
               Individual::mother_Individual,
               Individual::father_Individual

=cut

sub fetch_by_dbID {
  my $self = shift;
  my $dbID = shift;

  defined($dbID) || throw('dbID argument expected');

  my $sth = $self->prepare
    (q{SELECT i.individual_id, i.name, i.description,
              i.gender, i.father_individual_id, i.mother_individual_id
       FROM   individual i
       WHERE  i.individual_id = ?});

  $sth->execute($dbID);

  my $result = $self->_objs_from_sth($sth);

  $sth->finish();

  return undef if(!@$result);

  return $result->[0];
}



=head2 fetch_all_by_name

  Arg [1]    : string $name the name of the individuals to retrieve
  Example    : my @inds = @{$ind_adaptor->fetch_all_by_name('CEPH1332.05')};
  Description: Retrieves all individuals with a specified name.  Individual
               names may be non-unique which is why this method returns a
               reference to a list.
  Returntype : reference to a list of Individual ids
  Exceptions : throw if no argument passed
  Caller     : general

=cut

sub fetch_all_by_name {
  my $self = shift;
  my $name = shift;

  defined($name) || throw("name argument expected");

  my $sth = $self->prepare
    (q{SELECT i.individual_id, i.name, i.description,
              i.gender, i.father_individual_id, i.mother_individual_id
       FROM   individual i
       WHERE  i.name = ?});

  $sth->execute($name);

  my $result =  $self->_objs_from_sth($sth);

  $sth->finish();

  return $result;
}




=head2 fetch_all_by_Population

  Arg [1]    : Bio::EnsEMBL::Variation::Population $pop
  Example    : my $pop = $pop_adaptor->fetch_by_name('PACIFIC');
               foreach my $ind (@{$ia->fetch_all_by_Population($pop)}) {
                 print $ind->name(), "\n";
               }
  Description: Retrieves all individuals from a specified population 
  Returntype : reference to list of Bio::EnsEMBL::Variation::Individual objects
  Exceptions : throw if incorrect argument is passed
               warning if provided Population does not have an dbID
  Caller     : general

=cut

sub fetch_all_by_Population {
  my $self = shift;
  my $pop = shift;

  if(!ref($pop) || !$pop->isa('Bio::EnsEMBL::Variation::Population')) {
    throw("Bio::EnsEMBL::Variation::Population arg expected");
  }

  if(!$pop->dbID()) {
    warning("Population does not have dbID, cannot retrieve Individuals");
    return [];
  }

  my $sth = $self->prepare
    (q{SELECT i.individual_id, i.name, i.description,
              i.gender, i.father_individual_id, i.mother_individual_id
       FROM   individual i, individual_population ip
       WHERE  i.individual_id = ip.individual_id
       AND    ip.population_id = ?});

  $sth->execute($pop->dbID());

  my $results = $self->_objs_from_sth($sth);

  $sth->finish();

  return $results;
}




=head2 fetch_all_by_parent_Individual

  Arg [1]    : Bio::EnsEMBL::Variation::Individual
  Example    : my @children = @{$ia->fetch_all_by_parent_Individual($ind)};
  Description: Retrieves all individuals which are children of a provided
               parent individual.  This function operates under the assumptions
               that Male individuals can only be fathers, Female individuals
               can only be mothers and Unknown individuals can only be one
               or the other - not both.
  Returntype : reference to list of Bio::EnsEMBL::Variation::Individuals
  Exceptions : throw if incorrect argument passed
               warning if provided individual has no dbID 
  Caller     : general, Individual::get_all_child_Individuals

=cut

sub fetch_all_by_parent_Individual {
  my $self = shift;
  my $parent  = shift;

  if(!ref($parent) || !$parent->isa('Bio::EnsEMBL::Variation::Individual')) {
    throw("Bio::EnsEMBL::Variation::Individual argument expected");
  }

  if(!defined($parent->dbID())) {
    warning("Cannot fetch child Individuals for parent without dbID");
    return [];
  }

  my $gender = $parent->gender() || '';

  my $father_sql =
    q{SELECT i.individual_id, i.name, i.description,
             i.gender, i.father_individual_id, i.mother_individual_id
      FROM   individual i
      WHERE  i.father_individual_id = ?};
  my $mother_sql =
    q{SELECT i.individual_id, i.name, i.description,
              i.gender, i.father_individual_id, i.mother_individual_id
       FROM   individual i
       WHERE  i.mother_individual_id = ?};

  if($gender eq 'Male') {
    my $sth = $self->prepare($father_sql);
    $sth->execute($parent->dbID());
    my $result = $self->_objs_from_sth($sth);
    $sth->finish();
    return $result;
  }
  elsif($gender eq 'Female') {
    my $sth = $self->prepare($mother_sql);
    $sth->execute($parent->dbID());
    my $result = $self->_objs_from_sth($sth);
    $sth->finish();
    return $result;
  }

  # unknown gender

  my $sth = $self->prepare($mother_sql);
  $sth->execute($parent->dbID());
  my $result = $self->_objs_from_sth($sth);
  $sth->finish();

  # if this parent was a mother, finish now and return results
  return if(@$result);

  # otherwise assume was a father (or nothing)
  $sth = $self->prepare($father_sql);
  $sth->execute($parent->dbID());
  $result = $self->_objs_from_sth($sth);
  $sth->finish();

  return $result;
}





#
# private method, constructs Individuals from an executed statement handle
# ordering of columns must be consistant
#
sub _objs_from_sth {
  my $self = shift;
  my $sth = shift;

  my ($dbID, $name, $desc, $gender, $father_id, $mother_id);

  $sth->bind_columns(\$dbID, \$name, \$desc, \$gender,
                     \$father_id, \$mother_id);

  my %seen;
  my %wanted_fathers;
  my %wanted_mothers;

  my @inds;


  while($sth->fetch()) {
    # get objects for mother and father if they were already constructed
    # otherwise may have to be lazy-loaded later

    my $father;
    if(defined($father_id)) {
      $father = $seen{$father_id};
      if(!$father) {
        $wanted_fathers{$father_id} ||= [];
        push @{$wanted_fathers{$father_id}}, $father_id;
      }
    }
    my $mother;
    if(defined($mother_id)) {
      $mother = $seen{$mother_id};
      if(!$mother) {
        push @{$wanted_mothers{$mother_id}}, $mother_id;
      }
    }


    my $ind = $seen{$dbID} ||= Bio::EnsEMBL::Variation::Individual->new
      (-dbID        => $dbID,
       -adaptor     => $self,
       -description => $desc,
       -gender      => $gender,
       -name        => $name,
       -father_individual => $father,
       -mother_individual => $mother,
       -father_individual_id => $father_id,
       -mother_individual_id => $mother_id);

    $seen{$dbID} = $ind;

    push @inds, $ind;
  }

  # load any of the 'wanted' parent individuals that we did not have at the
  # of creation, but which we have now

  foreach my $wanted_id (keys %wanted_fathers) {
    if($seen{$wanted_id}) {
      # add father to every child that wanted it
      foreach my $ind_id (@{$wanted_fathers{$wanted_id}}) {
        $seen{$ind_id}->father_Individual($seen{$wanted_id});
      }
    }
  }
  foreach my $wanted_id (keys %wanted_mothers) {
    if($seen{$wanted_id}) {
      # add mother to every child that wanted it
      foreach my $ind_id (@{$wanted_mothers{$wanted_id}}) {
        $seen{$ind_id}->mother_Individual($seen{$wanted_id});
      }
    }
  }

  return \@inds;

}

sub _tables{return (['individual','i'])}

sub _columns{
    return qw(individual_id name description gender father_individual_id mother_individual_id);
}


1;
