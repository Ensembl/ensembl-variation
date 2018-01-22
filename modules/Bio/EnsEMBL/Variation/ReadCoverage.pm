=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut


=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <dev@ensembl.org>.

 Questions may also be sent to the Ensembl help desk at
 <helpdesk@ensembl.org>.

=cut

# Ensembl module for Bio::EnsEMBL::Variation::ReadCoverage
#
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
        -sample  => $sample);

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
    Bio::EnsEMBL::Variation::Sample - the sample in which the allele was recorded
    
  Example    :
    $rc = Bio::EnsEMBL::Variation::ReadCoverage->new
       (-start   => 100,
        -end     => 100,
        -slice   => $slice,
        -level  => 1,
        -sample => $sample);

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
  my ($level, $sample) =
    rearrange([qw(LEVEL SAMPLE)], @_);

  $self->{'level'}    = $level;
  $self->{'sample'}   = $sample;

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

  Arg [1]    : Bio::EnsEMBL::Variation::Sample $sample (optional)
               The new value to set the sample attribute to
  Example    : $sample = $rc->sample();
  Description: Getter/Setter for the sample attribute
  Returntype : Bio::EnsEMBL::Variation::Sample
  Exceptions : throw on incorrect argument
  Caller     : general
  Status     : Stable

=cut

sub sample {
  my $self = shift;
  my $sample = shift;
  if ($sample) {
    if(!ref($sample) || !$sample->isa('Bio::EnsEMBL::Variation::Sample')) {
      throw('Bio::EnsEMBL::Variation::Sample argument expected.');
    }
    $self->{'sample'} = $sample;
  }

  return $self->{'sample'};
}

1;
