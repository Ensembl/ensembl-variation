=head1 LICENSE

 Copyright (c) 1999-2012 The European Bioinformatics Institute and
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
        -sample  => $individual);

    $rc = $rc->transform('supercontig');

    print $rc->start(), "-", $rc->end(), "\n";


=head1 DESCRIPTION

This is a class representing the read coverage information
from the ensembl-variation database. A ReadCoverage behaves as any other Ensembl feature.

See B<Bio::EnsEMBL::Feature>.

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

  Arg [-SAMPLE] :
    Bio::EnsEMBL::Variation::Individual - the individual
                     in which the allele was recorded
    
  Example    :
    $rc = Bio::EnsEMBL::Variation::ReadCoverage->new
       (-start   => 100,
        -end     => 100,
        -slice   => $slice,
        -level  => 1,
        -sample => $individual);

  Description: Constructor. Instantiates a new ReadCoverage object.
  Returntype : Bio::EnsEMBL::Variation::ReadCoverage
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);
  my ($level, $individual) =
    rearrange([qw(LEVEL SAMPLE)], @_);

  $self->{'level'}    = $level;
  $self->{'sample'}   = $individual;

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
    Status      : Stable

=cut

sub level{
    my $self = shift;
    return $self->{'level'} = shift if (@_);
    return $self->{'level'};
}


=head2 sample

  Arg [1]    : Bio::EnsEMBL::Variation::Individual $newval (optional)
               The new value to set the sample attribute to
  Example    : $individual = $rc->sample();
  Description: Getter/Setter for the individual attribute
  Returntype : Bio::EnsEMBL::Variation::Individual
  Exceptions : throw on incorrect argument
  Caller     : general
  Status     : Stable

=cut

sub sample{
  my $self = shift;

  if(@_) {
    if(!ref($_[0]) || !$_[0]->isa('Bio::EnsEMBL::Variation::Individual')) {
      throw('Bio::EnsEMBL::Variation::Individual argument expected.');
    }
    $self->{'sample'} = shift;
  }

  return $self->{'sample'};
}

1;
