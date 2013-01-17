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

# Ensembl module for Bio::EnsEMBL::Variation::ReadCoverage
#
# Copyright (c) 2008 Ensembl
#


=head1 NAME

Bio::EnsEMBL::Variation::ReadCoverageCollection - A collection of coverage reagion for a read.

=head1 SYNOPSIS

    # Read coverage collection feature representing a genomic region covered by a read

    $rcc = Bio::EnsEMBL::Variation::ReadCoverageCollection->new
       (-start   => 100,
        -end     => 200,
        -strand  => 1,
        -slice   => $slice,
        -window_start => 1,
        -window_end => 1000,
        -window_size => 50,
        -read_coverage_avg => 100,
        -read_coverage_min => 0,
        -read_coverage_max => 200,
        -sample  => $individual);


    print $rcc->start(), "-", $rcc->end(), "\n";


=head1 DESCRIPTION

This is a class representing the read coverage collection information
from the ensembl-variation database. A ReadCoverageCollection behaves as any other Ensembl feature collection.
Object for storing read_coverage scores. The scores are averaged (also the minimum and maximum scores) over different window sizes to speed up drawing over large regions. The scores are packed as integer and stored in a string
See B<Bio::EnsEMBL::Feature>.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::ReadCoverageCollection;

use Bio::EnsEMBL::Feature;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument  qw(rearrange);

our @ISA = ('Bio::EnsEMBL::Feature');


=head2 new

  Arg [-ADAPTOR] :
    see superclass constructor
  Arg [-START] :
    start relative to the slice,see superclass constructor
  Arg [-END] :
    end relative to the slice,see superclass constructor
  Arg [-STRAND] :
    1 or -1, same as strand for the slice,see superclass constructor
  Arg [-SLICE] :
    see superclass constructor
  Arg [-SEQ_REGION_START] :
    start relative to chromosome,see superclass constructor
  Arg [-SEQ_REGION_END] :
    end relative to chromosome,see superclass constructor
  Arg [-SEQ_REGION_STRAND] :
    always 1, see superclass constructor
  int Arg [-SCORE_MIN] :
    average read_coverage for a window
  Arg [-SCORE_AVG] :
    minimum read_coverage for a window
  Arg [-SCORE_MAX] :
    maximum read_coverage for a window
  Arg [-SAMPLE_ID] :
    int - the individual in which the read_covarage is recorded
    
  Example    :
    $rc = Bio::EnsEMBL::Variation::ReadCoverage->new
       (-start   => 100,
        -end     => 100,
        -slice   => $slice,
        -window_size  => 50,
        -sample_id => 1);

  Description: Constructor. Instantiates a new ReadCoverageCollection object.
  Returntype : Bio::EnsEMBL::Variation::ReadCoverageCollection
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);
  my ($window_size,$window_start,$window_end,$sample_id,$read_coverage_avg,$read_coverage_min,$read_coverage_max,$y_axis_min,$y_axis_max) =
    rearrange([qw(WINDOW_SIZE WINDOW_START WINDOW_END SAMPLE_ID READ_COVERAGE_AVG READ_COVERAGE_MIN READ_COVERAGE_MAX Y_AXIS_MIN Y_AXIS_MAX)], @_);

  $self->{'window_size'}  = $window_size;
  $self->{'window_start'} = $window_start;
  $self->{'window_end'}   = $window_end;
  $self->{'sample_id'}   = $sample_id;
  $self->{'read_coverage_avg'} = $read_coverage_avg;
  $self->{'read_coverage_min'} = $read_coverage_min;
  $self->{'read_coverage_max'} = $read_coverage_max;
  $self->{'y_axis_min'} = $y_axis_min;
  $self->{'y_axis_max'} = $y_axis_max;
  return $self;
}

=head new_fast

  Arg [1]    : hash reference $hashref
  Example    : none
  Description: This is an ultra fast constructor which requires knowledge of
               the objects internals to be used.
  Returntype :
  Exceptions : none
  Caller     :

sub new_fast {

  my ($class, $hashref) = @_;

  return bless $hashref, $class;

}

=head2 window_size

    Arg[1]      : int $newval (optional)
                  The new value to set the window_size attribute to
    Example     : $window_size = $obj->window_size();
    Description : Getter/Setter for the window_size attribute.
                  the window size this feature has been seen in the genome
    ReturnType  : int 
    Exceptions  : none
    Caller      : general
    Status      : At Risk

=cut

sub window_size{
    my $self = shift;
    return $self->{'window_size'} = shift if (@_);
    return $self->{'window_size'};
}

=head2 window_start

    Arg[1]      : int $newval (optional)
                  The new value to set the window_start attribute to
    Example     : $window_start = $obj->window_start();
    Description : Getter/Setter for the window_start attribute.
                  the window start this feature has been seen in the genome
    ReturnType  : int 
    Exceptions  : none
    Caller      : general
    Status      : At Risk

=cut

sub window_start{
    my $self = shift;
    return $self->{'window_start'} = shift if (@_);
    return $self->{'window_start'};
}

=head2 window_end

    Arg[1]      : int $newval (optional)
                  The new value to set the window_end attribute to
    Example     : $depth = $obj->window_end();
    Description : Getter/Setter for the window_end attribute.
                  the window end this feature has been seen in the genome
    ReturnType  : int 
    Exceptions  : none
    Caller      : general
    Status      : At Risk

=cut

sub window_end{
    my $self = shift;
    return $self->{'window_end'} = shift if (@_);
    return $self->{'window_end'};
}

=head2 read_coverage_avg

    Arg[1]      : int $newval (optional)
                  The new value to set the read_coverage_avg attribute to
    Example     : $avg = $obj->read_coverage_avg();
    Description : Getter/Setter for the score_avg attribute.
                  the average read_coverage this feature has
    ReturnType  : int 
    Exceptions  : none
    Caller      : general
    Status      : At Risk

=cut

sub read_coverage_avg{
    my $self = shift;
    return $self->{'read_coverage_avg'} = shift if (@_);
    return $self->{'read_coverage_avg'};
}

=head2 read_coverage_min

    Arg[1]      : int $newval (optional)
                  The new value to set the read_coverage_min attribute to
    Example     : $min = $obj->read_coverage_min();
    Description : Getter/Setter for the score_min attribute.
                  the minimum read_coverage this feature has
    ReturnType  : int 
    Exceptions  : none
    Caller      : general
    Status      : At Risk

=cut

sub read_coverage_min{
    my $self = shift;
    return $self->{'read_coverage_min'} = shift if (@_);
    return $self->{'read_coverage_min'};
}

=head2 read_coverage_max

    Arg[1]      : int $newval (optional)
                  The new value to set the read_coverage_max attribute to
    Example     : $max = $obj->read_coverage_max();
    Description : Getter/Setter for the read_coverage_max attribute.
                  the maximum read_coverage this feature has
    ReturnType  : int 
    Exceptions  : none
    Caller      : general
    Status      : At Risk

=cut

sub read_coverage_max{
    my $self = shift;
    return $self->{'read_coverage_max'} = shift if (@_);
    return $self->{'read_coverage_max'};
}


=head2 sample_id

    Arg[1]      : int $newval (optional)
                  The new value to set the sample_id attribute to
    Example     : $sample_id = $obj->sample_id();
    Description : Getter/Setter for the individual dbId attribute.
                  the individual dbId this feature has been seen in the genome
    ReturnType  : int 
    Exceptions  : none
    Caller      : general
    Status      : At Risk

=cut

sub sample_id{
    my $self = shift;
    return $self->{'sample_id'} = shift if (@_);
    return $self->{'sample_id'};
}

=head2 y_axis_min

    Arg[1]      : int $newval (optional)
                  The new value to set the y_axiss_min attribute to
    Example     : $y_axis_min = $obj->y_axis_min();
    Description : Getter/Setter for the minimum of read_coverage for the collection of feature.

    ReturnType  : int 
    Exceptions  : none
    Caller      : general
    Status      : At Risk

=cut

sub y_axis_min{
    my $self = shift;
    return $self->{'y_axis_min'} = shift if (@_);
    return $self->{'y_axis_min'};
}

=head2 y_axis_max

    Arg[1]      : int $newval (optional)
                  The new value to set the y_axiss_max attribute to
    Example     : $y_axis_max = $obj->y_axis_max();
    Description : Getter/Setter for the maximum of read_coverage for the collection of feature.

    ReturnType  : int 
    Exceptions  : none
    Caller      : general
    Status      : At Risk

=cut

sub y_axis_max{
    my $self = shift;
    return $self->{'y_axis_max'} = shift if (@_);
    return $self->{'y_axis_max'};
}


1;
