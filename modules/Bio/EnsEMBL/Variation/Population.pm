# Ensembl module for Bio::EnsEMBL::Variation::Population
#
# Copyright (c) 2004 Ensembl
#


=head1 NAME

Bio::EnsEMBL::Variation::Population - A population represents a phenotypic group, ethnic group, set of individuals used in an assay, etc.

=head1 SYNOPSIS

    # Population
    $pop = Bio::EnsEMBL::Variation::Population->new
       (-name => 'WEST AFRICA',
        -description => 'Sub-Saharan Nations bordering Atlantic north' .
                        ' of Congo River, and Central/Southern Atlantic' .
                        ' Island Nations.');

    ...

    # print out all sub populations of a population
    # same could work for super populations

    print_sub_pops($pop);

    sub print_sub_pops {
      my $pop = shift;
      my $level = shift || 0;

      my $sub_pops = $pop->get_all_sub_Populations();

      foreach my $sp (@$sub_pops) {
        print ' ' x $level++,
              'name: ', $sp->name(),
              'desc: ', $sp->description(),
              'size: ', $sp->size(),"\n";
        print_sub_pops($sp, $level);
      }
    }



=head1 DESCRIPTION

This is a class representing a population.  A population may consist of any
grouping of individuals, including phenotypic groups (e.g. people with
diabetes), ethnic groups (e.g. caucasians), individuals used in an assay
(e.g. subjects in experiment X), etc.

Populations may be arranged into an arbitrary hierarchy of sub and super
populations.

=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::Population;

use Bio::EnsEMBL::Storable;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);

our @ISA = ('Bio::EnsEMBL::Storable');


=head2 new

  Arg [-dbID]: int - unique internal identifier of the population
  Arg [-ADAPTOR]: Bio::EnsEMBL::PopulationAdaptor
  Arg [-NAME]: string - name of the population
  Arg [-DESCRIPTION]: string - description of the population
  Arg [-SIZE]: int - the size of the population
  Arg [-SUB_POPULATIONS]: listref of Bio::EnsEMBL::Population objects 
  Example    : $pop = Bio::EnsEMBL::Variation::Population->new
       (-name => 'WEST AFRICA',
        -description => 'Sub-Saharan Nations bordering Atlantic north' .
                        ' of Congo River, and Central/Southern Atlantic' .
                        ' Island Nations.'
        -sub_populations => \@sub_pops);
  Description: Constructor. Instantiates a new Population object
  Returntype : Bio::EnsEMBL::Variation::Population
  Exceptions : none
  Caller     : general

=cut

sub new {
  my $caller = shift;

  my $class = ref($caller) || $caller;

  my ($dbID, $adaptor, $name, $desc, $size, $sub_pops) =
    rearrange(['DBID','ADAPTOR','NAME', 'DESCRIPTION', 'SIZE',
               'SUB_POPULATIONS'], @_);

  return bless {'dbID'        => $dbID,
                'adaptor'     => $adaptor,
                'name'        => $name,
                'description' => $desc,
                'size'        => $size,
                'sub_populations' => $sub_pops || []}, $class;
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



=head2 size

  Arg [1]    : int $newval (optional) 
               The new value to set the size attribute to
  Example    : $size = $obj->size()
  Description: Getter/Setter for the size attribute
  Returntype : int
  Exceptions : none
  Caller     : general

=cut

sub size{
  my $self = shift;
  return $self->{'size'} = shift if(@_);
  return $self->{'size'};
}



=head2 get_all_sub_Populations

  Arg [1]    : none
  Example    : foreach my $sub_pop (@{$pop->get_all_sub_Populations}) {
                 my $sub_sub_pops = $sub_pop->get_all_sub_Populations();
               }
  Description: Retrieves all populations which are conceptually a sub set
               of this population.
  Returntype : reference to list of Bio::EnsEMBL::Variation::Population objects
  Exceptions : none
  Caller     : general

=cut

sub get_all_sub_Populations {
  my $self = shift;

  ### TBD - lazy-load from database if not set
  return $self->{'sub_populations'};
}



=head2 get_all_super_Populations

  Arg [1]    : none
  Example    : foreach my $sup_pop (@{$pop->get_all_super_Populations}) {
                 my $sup_sup_pops = $sup_pop->get_all_super_Populations();
               }
  Description: Retrieves all populations which this population is a part of
               from the database.
               Super populations may not be directly added in order to avoid
               circular references and memory leaks.  You must add
               sub_Populations instead and store this in the database.
  Returntype : reference to list of Bio::EnsEMBL::Variation::Population objects
  Exceptions : none
  Caller     : general

=cut

sub get_all_super_Populatons {
  my $self = shift;

  ### TBD - lazy load from database - do not set to avoid circular references!
}



=head2 add_sub_Population

  Arg [1]    : Bio::EnsEMBL::Variation::Population $pop
  Example    : $pop->add_sub_Population($sub_pop);
               $sub_pop->add_super_Population($pop);
  Description: Adds a sub population to this population.
  Returntype : none
  Exceptions : throw on incorrect argument
  Caller     : general

=cut

sub add_sub_Population {
  my $self = shift;
  my $pop = shift;

  if(!ref($pop) || !$pop->isa('Bio::EnsEMBL::Variation::Population')) {
    throw('Bio::EnsEMBL::Variation::Population argument expected.');
  }

  if($pop == $self) {
    throw("Cannot add self as sub population.");
  }

  push @{$self->{'sub_populations'}}, $pop;

  return $pop;
}


1;
