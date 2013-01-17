=head1 LICENSE

 Copyright (c) 1999-2013 The European Bioinformatics Institute and
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
of one population or several.  A pedigree may be constructed using the father_Individual
and mother_Individual attributes.

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::Variation::Individual;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Variation::Sample;

our @ISA = ('Bio::EnsEMBL::Variation::Sample');

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
  Arg [-FATHER_INDIVIDUAL] :
    Bio::EnsEMBL::Variation::Individual - the father of this individual
  Arg [-MOTHER_INDIVIDUAL] :
    Bio::EnsEMBL::Variation::Individual - the mother of this individual
  Arg [-MOTHER_INDIVIDUAL_SAMPLE_ID] :
    int - set the internal id of the mother individual so that the actual
    mother Individual object can be retrieved on demand.
  Arg [-FATHER_INDIVIDUAL_SAMPLE_ID]:
    int - set the internal id of the mother individual so that the actual
    mother Individual object can be retrieved on demand.
  Arg [-TYPE_INDIVIDUAL]:
    int - name for the type of the individual (fully or partly inbred, outbred or mutant
  Arg [-TYPE_DESCRIPTION]:
    string - description of the type of individual
  Example    : $individual = Bio::EnsEMBL::Variation::Individual->new
                 (-name => 'WI530.07',
                  -description => 'african',
                  -gender => 'Male',
                  -father_individual => $father_ind,
                  -mother_individual => $mother_ind,
		  -type_individual => 'outbred',
		  -type_description => 'a single organism which breeds freely');
  Description: Constructor Instantiates an Individual object.
  Returntype : Bio::EnsEMBL::Variation::Individual
  Exceptions : throw if gender arg is provided but not valid
  Caller     : general
  Status     : Stable

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my ($dbID, $adaptor, $name, $desc, $display_flag, $gender, $father, $mother, $type_name, $type_desc,
      $father_id, $mother_id) =
    rearrange([qw(dbID adaptor name description display gender
                  father_individual mother_individual
		  type_individual type_description
                  father_individual_sample_id mother_individual_sample_id)], @_);

  if(defined($gender)) {
    $gender = ucfirst(lc($gender));
    if($gender ne 'Male' && $gender ne 'Female' && $gender ne 'Unknown') {
      throw('Gender must be one of "Male","Female","Unknown"');
    }
  }

  if (defined($type_name)){
      $type_name = ucfirst(lc($type_name));
      if ($type_name ne 'Fully_inbred' && $type_name ne 'Partly_inbred' && $type_name ne 'Outbred' && $type_name ne 'Mutant'){
	  throw('Type of individual must of one of: "fully_inbred", "partly_inbred", "outbred", "mutant"');
      }
  }
  
  return bless {
	'dbID'    => $dbID,
	'adaptor' => $adaptor,
	'name'    => $name,
	'description' => $desc,
	'display' => $display_flag,
	'gender'  => $gender,
	'father_individual' => $father,
	'mother_individual' => $mother,
	'type_individual' => $type_name,
	'type_description' => $type_desc,
	'_mother_individual_sample_id' => $mother_id,
	'_father_individual_sample_id' => $father_id,
  }, $class;
}


=head2 type_individual

    Arg [1]     : int $newval (optional)
                  The new value to set the type_individual attribute to                  
    Example     : $type_individual = $obj->type_individual();
    Description : Getter/Setter for the type_individual attribute
    Returntype  : int
    Exceptions  : none
    Caller      : general
    Status      : At Risk

=cut

sub type_individual{
    my $self = shift;
    if (@_){
	my $new_name = shift;
	return $self->{'type_individual'} = $new_name;
    }
    return $self->{'type_individual'};
}

=head2 type_description

    Arg [1]     : int $newval (optional)
                  The new value to set the type_description attribute to                  
    Example     : $type_description = $obj->type_description();
    Description : Getter/Setter for the type_description attribute
    Returntype  : int
    Exceptions  : none
    Caller      : general
    Status      : At Risk

=cut

sub type_description{
    my $self = shift;
    if (@_){
	my $new_desc = shift;
	return $self->{'type_description'} = $new_desc;
    }
    return $self->{'type_description'};
}

=head2 gender

  Arg [1]    : string $newval (optional)
               The new value to set the gender attribute to
  Example    : $gender = $obj->gender()
  Description: Getter/Setter for the gender attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

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


=head2 display

  Arg [1]    : string $newval (optional)
               The new value to set the display attribute to
  Example    : $display = $obj->display()
  Description: Getter/Setter for the display attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub display{
  my $self = shift;
  
  if(@_) {
    my $display = uc(shift);
    
    if($display ne 'UNDISPLAYABLE' && $display ne 'REFERENCE' && $display ne 'DISPLAYABLE' && $display ne 'DEFAULT') {
      throw('Display flag must be one of "REFERENCE", "DEFAULT", "DISPLAYABLE", "UNDISPLAYABLE"');
    }
    
    $self->{'display'} = $display;
  }
  
  return $self->{'display'};
}


=head2 get_all_Populations

   Args        : none
   Example     : $pops = $ind->get_all_Populations();
   Description : Getter for the Populations for this Individual. Returns
                 empty list if there are none.
   ReturnType  : listref of Bio::EnsEMBL::Population
   Exceptions  : none
   Caller      : general
   Status      : Stable

=cut

sub get_all_Populations{
  my $self = shift;
  
  if(!defined($self->{populations})) {
	if (defined ($self->{'adaptor'})){
	  my $pop_adaptor = $self->{'adaptor'}->db()->get_PopulationAdaptor();
	  $self->{populations} = $pop_adaptor->fetch_all_by_Individual($self);
	}
  }
  
  return $self->{populations};
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
  Status     : Stable

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
     defined($self->{'_father_individual_sample_id'})) {
    $self->{'father_individual'} =
      $self->adaptor->fetch_by_dbID($self->{'_father_individual_sample_id'});
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
  Status     : Stable

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
     defined($self->{'_mother_individual_sample_id'})) {
    $self->{'mother_individual'} =
      $self->adaptor->fetch_by_dbID($self->{'_mother_individual_sample_id'});
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
  Status     : Stable

=cut

sub get_all_child_Individuals {
  my $self = shift;

  if(!$self->adaptor()) {
    warning("Cannot retrieve child individuals without attached adaptor.");
  }

  return $self->adaptor()->fetch_all_by_parent_Individual($self);
  
}

1;
