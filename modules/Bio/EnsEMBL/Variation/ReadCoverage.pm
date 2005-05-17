# Ensembl module for Bio::EnsEMBL::Variation::ReadCoverage
#
# Copyright (c) 2005 Ensembl
#


=head1 NAME

Bio::EnsEMBL::Variation::ReadCoverage - A coverage reagion for a read.

=head1 SYNOPSIS

    # Read coverage feature representing a genomic region covered by 1 read

    $rc = Bio::EnsEMBL::Variation::ReadCoverage->new
       (-start   => 100,
        -end     => 200,
        -slice   => $slice,
        -level   => 1.
        -population => $population);

    $rc = $rc->transform('supercontig');

    print $rc->start(), "-", $rc->end(), "\n";


=head1 DESCRIPTION

This is a class representing the read coverage information
from the ensembl-variation database. A ReadCoverage behaves as any other Ensembl feature.

See B<Bio::EnsEMBL::Feature>.

=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::ReadCoverage;

use Bio::EnsEMBL::Feature;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument  qw(rearrange);

our @ISA = ('Bio::EnsEMBL::Feature');


=head2 new

  Arg [-ADAPTOR] :
    see superclass constructor

  Arg [-START] :
    see superclass constructor
  Arg [-END] :
    see superclass constructor
  Arg [-SLICE] :
    see superclass constructor

  Arg [-LEVEL] :
    int - the number of times the region represented by start and end has been seen

  Arg [-POPULATION] :
    Bio::EnsEMBL::Variation::Population - the population
                     in which the allele was recorded
    
  Example    :
    $rc = Bio::EnsEMBL::Variation::ReadCoverage->new
       (-start   => 100,
        -end     => 100,
        -slice   => $slice,
        -level  => 1,
        -population => $population);

  Description: Constructor. Instantiates a new ReadCoverage object.
  Returntype : Bio::EnsEMBL::Variation::ReadCoverage
  Exceptions : none
  Caller     : general

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);
  my ($level, $population) =
    rearrange([qw(LEVEL POPULATION)], @_);

  $self->{'level'}    = $level;
  $self->{'population'}   = $population;

  return $self;
}


=head2 level

    Arg[1]      : int $newval (optional)
                  The new value to set the level attribute to
    Example     : $depth = $obj->level();
    Description : Getter/Setter for the level attribute. The level is
                  the number of times this feature has been seen in the genome
    ReturnType  : int 
    Exceptions  : none
    Caller      : general

=cut

sub level{
    my $self = shift;
    return $self->{'level'} = shift if (@_);
    return $self->{'level'};
}


=head2 population

  Arg [1]    : Bio::EnsEMBL::Variation::Population $newval (optional)
               The new value to set the population attribute to
  Example    : $population = $rc->population();
  Description: Getter/Setter for the population attribute
  Returntype : Bio::EnsEMBL::Variation::Population
  Exceptions : throw on incorrect argument
  Caller     : general

=cut

sub population{
  my $self = shift;

  if(@_) {
    if(!ref($_[0]) || !$_[0]->isa('Bio::EnsEMBL::Variation::Population')) {
      throw('Bio::EnsEMBL::Variation::Population argument expected.');
    }
    $self->{'population'} = shift;
  }

  return $self->{'population'};
}

1;
