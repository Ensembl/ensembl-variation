# Ensembl module for Bio::EnsEMBL::Variation::Individual
#
# Copyright (c) 2004 Ensembl
#


=head1 NAME

Bio::EnsEMBL::Variation::Individual - A single member of a population.

=head1 SYNOPSIS

    $individual = Bio::EnsEMBL::Variation::Individual->new
       (-name => 'WI530.07',
        -description => 'african',
        -gender => 'Male',
        -population => $nigerian_pop,
        -father_individual => $father_ind,
        -mother_individual => $mother_ind);

    ...

    print $individual->name(), ' - ', $individual->description(), "\n";
    print "Gender: ", $individual->gender(), "\n";
    print $individual->mother_Individual->name()
       if($individual->mother_Individual());
    print $individual->father_Individual->name()
       if($individual->father_Individual());



=head1 DESCRIPTION

This is a class representing a single individual.  An individual may be part
of a population.  A pedigree may be constructed using the father_Individual
and mother_Individual attributes.

=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::Variation::Individual;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Storable;

our @ISA = ('Bio::EnsEMBL::Storable');

=head2 new

  Arg [-dbID] :
    int - unique internal identifier
  Arg [-ADAPTOR] :
    Bio::EnsEMBL::Variation::DBSQL::IndividualAdaptor
  Arg [-NAME] :
    string - name of this individual
  Arg [-DESCRIPTION] :
    string description - description of this individual
  Arg [-GENDER] :
    string - must be one of 'Male', 'Female', 'Unknown'
  Arg [-POPULATION] :
    Bio::EnsEMBL::Variation::Population - the pop this individual is from
  Arg [-FATHER_INDIVIDUAL] :
    Bio::EnsEMBL::Variation::Individual - the father of this individual
  Arg [-MOTHER_INDIVIDUAL] :
    Bio::EnsEMBL::Variation::Individual - the mother of this individual
  Arg [-MOTHER_INDIVIDUAL_ID] :
    int - set the internal id of the mother individual so that the actual
    mother Individual object can be retrieved on demand.
  Arg [-FATHER_INDIVIDUAL_ID]:
    int - set the internal id of the mother individual so that the actual
    mother Individual object can be retrieved on demand.
  Example    : $individual = Bio::EnsEMBL::Variation::Individual->new
                 (-name => 'WI530.07',
                  -description => 'african',
                  -gender => 'Male',
                  -population => $nigerian_pop,
                  -father_individual => $father_ind,
                  -mother_individual => $mother_ind);
  Description: Constructor Instantiates an Individual object.
  Returntype : Bio::EnsEMBL::Variation::Individual
  Exceptions : throw if gender arg is provided but not valid
  Caller     : general

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my ($dbID, $adaptor, $name, $desc, $gender, $pop, $father, $mother,
      $father_id, $mother_id) =
    rearrange([qw(dbID adaptor name description gender population
                  father_individual mother_individual
                  father_individual_id mother_individual_id)], @_);

  if(defined($gender)) {
    $gender = ucfirst(lc($gender));
    if($gender ne 'Male' && $gender ne 'Female' && $gender ne 'Unknown') {
      throw('Gender must be one of "Male","Female","Unknown"');
    }
  }

  return bless {'dbID'    => $dbID,
                'adaptor' => $adaptor,
                'name'    => $name,
                'description' => $desc,
                'gender'  => $gender,
                'population' => $pop,
                'father_individual' => $father,
                'mother_individual' => $mother,
                '_mother_individual_id' => $mother_id,
                '_father_individual_id' => $father_id}, $class;
}



=head2 name

  Arg [1]    : string $newval (optional)
               The new value to set the name attribute to
  Example    : $name = $obj->name()
  Description: Getter/Setter for the name attribute
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub name{
  my $self = shift;
  return $self->{'name'} = shift if(@_);
  return $self->{'name'};
}



=head2 description

  Arg [1]    : string $newval (optional)
               The new value to set the description attribute to
  Example    : $description = $obj->description()
  Description: Getter/Setter for the description attribute
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub description{
  my $self = shift;
  return $self->{'description'} = shift if(@_);
  return $self->{'description'};
}



=head2 gender

  Arg [1]    : string $newval (optional)
               The new value to set the gender attribute to
  Example    : $gender = $obj->gender()
  Description: Getter/Setter for the gender attribute
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub gender{
  my $self = shift;

  if(@_) {
    my $gender = ucfirst(lc(shift));

    if($gender ne 'Male' && $gender ne 'Female' && $gender ne 'Unknown') {
      throw('Gender must be one of "Male","Female","Unknown"');
    }
    $self->{'gender'} = $gender;
  }

  return $self->{'gender'};
}



=head2 population

  Arg [1]    : string $newval (optional) 
               The new value to set the population attribute to
  Example    : $population = $obj->population()
  Description: Getter/Setter for the population attribute
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub population{
  my $self = shift;

  if(@_) {
    my $pop = shift;
    if(defined($pop) && (!ref($pop) ||
                         !$pop->isa('Bio::EnsEMBL::Variation::Population'))) {
      throw('Bio::EnsEMBL::Variation::Population arg expected');
    }
    $self->{'population'} = $pop;
  }

  return $self->{'population'};
}



=head2 father_Individual

  Arg [1]    : string $newval (optional)
               The new value to set the father_Individual attribute to
  Example    : $father_Individual = $obj->father_Individual()
  Description: Getter/Setter for the father of this Individual. If this
               has not been set manually and this Individual has an attached
               adaptor, an attempt will be made to lazy-load it from the
               database.
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub father_Individual{
  my $self = shift;

  if(@_) {
    my $ind = shift;
    if(defined($ind) && (!ref($ind) ||
                         !$ind->isa('Bio::EnsEMBL::Variation::Individual'))) {
      throw('Bio::EnsEMBL::Variation::Individual arg expected');
    }
    if($ind->gender() eq 'Female') {
      throw("Father individual may not have gender of Female");
    }
    return $self->{'father_individual'} = $ind;
  }

  # lazy-load mother if we can
  if(!defined($self->{'father_individual'}) && $self->adaptor() &&
     defined($self->{'_father_individual_id'})) {
    $self->{'father_individual'} =
      $self->adaptor->fetch_by_dbID($self->{'_father_individual_id'});
  }

  return $self->{'father_individual'};
}



=head2 mother_Individual

  Arg [1]    : string $newval (optional) 
               The new value to set the mother_Individual attribute to
  Example    : $mother_Individual = $obj->mother_Individual()
  Description: Getter/Setter for the mother of this individual. If this
               has not been set manually and this Individual has an attached
               adaptor, an attempt will be made to lazy-load it from the
               database.
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub mother_Individual{
  my $self = shift;

  if(@_) {
    my $ind = shift;
    if(defined($ind) && (!ref($ind) ||
                         !$ind->isa('Bio::EnsEMBL::Variation::Individual'))) {
      throw('Bio::EnsEMBL::Variation::Individual arg expected');
    }
    if($ind->gender() eq 'Male') {
      throw("Mother individual may not have gender of Male");
    }
    return $self->{'mother_individual'} = $ind;
  }

  # lazy-load mother if we can
  if(!defined($self->{'mother_individual'}) && $self->adaptor() &&
     defined($self->{'_mother_individual_id'})) {
    $self->{'mother_individual'} =
      $self->adaptor->fetch_by_dbID($self->{'_mother_individual_id'});
  }

  return $self->{'mother_individual'};
}



=head2 get_all_child_Individuals

  Arg [1]    : none
  Example    : foreach my $c (@{$ind->get_all_child_Individuals}) {
                 print "Child: " $c->name(), "\n";
               }
  Description: Retrieves all individuals from the database which are children
               of this individual. This will only work if this Individual
               object has been stored in the database and has an attached
               adaptor.
  Returntype : reference to list of Bio::EnsEMBL::Variation::Individual objects
  Exceptions : warning if this object does not have an attached adaptor
  Caller     : general

=cut

sub get_all_child_Individuals {
  my $self = shift;

  if(!$self->adaptor()) {
    warning("Cannot retrieve child individuals without attached adaptor.");
  }

  ### TODO: implement this.  DO NOT CACHE because would create circ. ref and
  ### memory leak.
}



1;
