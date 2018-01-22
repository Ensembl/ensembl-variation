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
 developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

 Questions may also be sent to the Ensembl help desk at
 <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Variation::TranscriptDiplotype

=head1 SYNOPSIS

  use Bio::EnsEMBL::Variation::TranscriptDiplotype;

=head1 DESCRIPTION

Class to represent an individual/sample's haplotype
pair across a transcript.

=cut

package Bio::EnsEMBL::Variation::TranscriptDiplotype;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Digest::MD5 qw(md5_hex);

=head2 new

  Arg [-haplotypes] : arrayref of Bio::EnsEMBL::Variation::CDSHaplotypes
  Arg [-samples]    : arrayref of sample names

  Example           : my $ch = Bio::EnsEMBL::Variation::TranscriptDiplotype->new(
                        -HAPLOTYPES => [$ht1, $ht2],
                        -SAMPLES    => [$sample_name1, $sample_name2],
                      );

  Description: Constructor.  Instantiates a new TranscriptDiplotype object.
  Returntype : Bio::EnsEMBL::Variation::DBSQL::TranscriptDiplotype
  Exceptions : none
  Caller     : Bio::EnsEMBL::Variation::TranscriptHaplotypeContainer
  Status     : at risk

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  my ($haplotypes, $samples) = rearrange([qw(HAPLOTYPES SAMPLES)], @_);

  assert_ref($haplotypes, 'ARRAY', 'Haplotypes listref');
  assert_ref($_, 'Bio::EnsEMBL::Variation::TranscriptHaplotype', 'First member of haplotypes listref') for @$haplotypes;
  
  my $self = {
    _haplotypes => $haplotypes,
    samples => $samples,
  };

  bless($self, $class);
  
  return $self;
}


=head2 name

  Example    : my $name = $td->name()
  Description: Get a name for this diplotype based on the transcript identifier
               and a pair of concatenated strings consisting of the constituent
               haplotypes' differences to the reference
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub name {
  my $self = shift;
  
  if(!exists($self->{name})) {
    my @names = map {$_->name()} @{$self->haplotypes()};

    my $ref_part = (split(':', $names[0]))[0];

    $self->{name} = $ref_part.':'.join("_", map {(split(':', $_))[1]} @names);
  }
  
  return $self->{name};
}


=head2 haplotypes

  Example    : my $haplotypes = $td->haplotypes()
  Description: Gets the constituent haplotypes of this diplotype
  Returntype : arrayref of Bio::EnsEMBL::Variation::TranscriptHaplotypes
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub haplotypes {
  return $_[0]->{_haplotypes};
}


=head2 type

  Example    : my $type = $td->type()
  Description: Get the type of the constituent haplotypes of this TranscriptDiplotype,
               one of "CDS" or "Protein"
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub type {
  my $self = shift;
  return $self->{type} ||= $self->haplotypes->[0]->type;
}


=head2 container

  Example    : my $container = $td->container()
  Description: Get the TranscriptHaplotypeContainer object used to generate this
               TranscriptDiplotype
  Returntype : Bio::EnsEMBL::TranscriptHaplotypeContainer
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub container {
  return $_[0]->haplotypes->[0]->container;
}


=head2 transcript

  Example    : my $tr = $td->transcript()
  Description: Get the Transcript object used as the reference for this
               TranscriptDiplotype
  Returntype : Bio::EnsEMBL::Transcript
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub transcript {
  return $_[0]->haplotypes->[0]->transcript;
}

=head2 count

  Example    : my $count = $td->count()
  Description: Counts the number of times this diplotype has been observed
               within the individuals defined in the container
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub count {
  my $self = shift;
  
  if(!exists($self->{count})) {
    $self->{count} = scalar @{$self->{samples}};
  }
  
  return $self->{count};
}


=head2 frequency

  Example    : my $frequency = $td->frequency()
  Description: Get the observed frequency of this diplotype amongst
               those defined in the container
  Returntype : float
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub frequency {
  my $self = shift;
  
  if(!exists($self->{frequency})) {
    my $total = $self->container->_total_sample_count;
    my $count = $self->count;
    $self->{frequency} = $total > 0 ? $count / $total : 0;
  }
  
  return $self->{frequency};
}


=head2 get_all_population_counts

  Example    : my %counts = %{$td->get_all_population_counts()}
  Description: Counts the number of times this diplotype has been observed
               in each population defined in the container
  Returntype : hashref e.g. { $population_name => $count }
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_population_counts {
  my $self = shift;
  
  if(!exists($self->{population_counts})) {
    my $sample_pop_hash = $self->container->_get_sample_population_hash();
    my $counts = {};
    
    foreach my $sample(@{$self->{samples}}) {
      $counts->{$_}++ for keys %{$sample_pop_hash->{$sample}};
    }
    
    $self->{population_counts} = $counts;
  }
  
  return $self->{population_counts};
}


=head2 get_all_population_frequencies

  Example    : my %freqs = %{$td->get_all_population_frequencies()}
  Description: Gets the observed frequency of this diplotype in each population
               defined in the container
  Returntype : hashref e.g. { $population_name => $freq }
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_population_frequencies {
  my $self = shift;
  
  if(!exists($self->{population_frequencies})) {
    my $totals = $self->container->_sample_counts_by_population;
    my $counts = $self->get_all_population_counts;
    
    my %freqs = map {$_ => $counts->{$_} / $totals->{$_}} keys %$counts;
    
    $self->{population_frequencies} = \%freqs;
  }
  
  return $self->{population_frequencies};
}

sub _hex {
  my $self = shift;
  return $self->{hex} ||= md5_hex(join("_", map {$_->_hex} @{$self->haplotypes}));
}

## Convert this object to a hash that can be written as JSON.
## Basically just deletes "private" keys starting with "_"
sub TO_JSON {
  my $self = shift;

  $self->$_ for qw(type name count frequency get_all_population_frequencies);
  
  # make a hash copy of self
  my %copy = %{$self};

  # copy haplotypes as hexes
  $copy{haplotypes} = [map {$_->_hex} @{$self->haplotypes}];
  
  # delete keys starting with _
  delete $copy{$_} for grep {$_ =~ /^\_/} keys %copy;
  
  return \%copy;
}

1;
