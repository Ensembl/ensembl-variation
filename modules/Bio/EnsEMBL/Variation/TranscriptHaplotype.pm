=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Variation::TranscriptHaplotype

=head1 SYNOPSIS

  use Bio::EnsEMBL::Variation::TranscriptHaplotype;

=head1 DESCRIPTION

A helper class for representing a transcript sequence modified by sample
genotypes. Not to be used directly.

=cut

package Bio::EnsEMBL::Variation::TranscriptHaplotype;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Argument qw(rearrange);

use Bio::AlignIO;
use Bio::SimpleAlign;
use Bio::LocatableSeq;

=head2 new

  Arg [-CONTAINER]:  Bio::EnsEMBL::Variation::TranscriptHaplotypeContainer
  Arg [-TYPE]:       string - "cds" or "protein"
  Arg [-SEQ]:        string
  Arg [-HEX]:        string
  Arg [-INDEL]:      bool

  Example    : my $ch = Bio::EnsEMBL::Variation::TranscriptHaplotype->new(
                  -CONTAINER => $container,
                  -TYPE      => 'protein',
                  -SEQ       => $seq,
                  -HEX       => $hex,
                  -INDEL     => $indel
               );

  Description: Constructor.  Instantiates a new TranscriptHaplotype object.
  Returntype : Bio::EnsEMBL::Variation::DBSQL::TranscriptHaplotype
  Exceptions : none
  Caller     : internal
  Status     : Stable

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  my ($hex, $type, $seq, $indel, $container) = rearrange([qw(HEX TYPE SEQ INDEL CONTAINER)], @_);
  
  my $self = {
    type => $type,
    seq => $seq,
    has_indel => $indel,
    hex => $hex,
    other_hexes => {},
    _container => $container
  };
  
  bless($self, $class);
  
  return $self;
}


=head2 name

  Example    : my $name = $th->name()
  Description: Get a name for this haplotype based on the transcript identifier
               and a concatenated string of its difference to the reference
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub name {
  my $self = shift;
  
  if(!exists($self->{name})) {
    $self->{name} = sprintf("%s:%s", $self->_reference_name, join(",", @{$self->_get_raw_diffs}) || 'REF');
  }
  
  return $self->{name};
}


=head2 container

  Arg [1]    : Bio::EnsEMBL::Variation::TranscriptHaplotypeContainer $container (optional)
               Set the container for this TranscriptHaplotype
  Example    : my $container = $th->container()
  Description: Getter/Setter for the container object.
  Returntype : Bio::EnsEMBL::Variation::TranscriptHaplotypeContainer
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub container {
  my $self = shift;
  $self->{_container} = shift if @_;
  return $self->{_container};
}


=head2 type

  Arg [1]    : string $type (optional)
  Example    : my $type = $th->type()
  Description: Getter/Setter for the type of this TranscriptHaplotype,
               one of "CDS" or "Protein"
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub type {
  my $self = shift;
  $self->{type} = shift if @_;
  return $self->{type};
}


=head2 seq

  Arg [1]    : string $seq (optional)
  Example    : my $seq = $th->seq()
  Description: Getter/Setter for the sequence of this TranscriptHaplotype
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub seq {
  my $self = shift;
  $self->{seq} = shift if @_;
  return $self->{seq};
}


=head2 has_indel

  Arg [1]    : bool $has_indel (optional)
  Example    : my $has_indel = $th->has_indel()
  Description: Getter/Setter for the flag indicating if this
               TranscriptHaplotype has an insertion or deletion
               relative to the reference
  Returntype : bool
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub has_indel {
  my $self = shift;
  $self->{has_indel} = shift if @_;
  return $self->{has_indel};
}


=head2 transcript

  Example    : my $tr = $th->transcript()
  Description: Get the Transcript object used as the reference for this
               TranscriptHaplotype
  Returntype : Bio::EnsEMBL::Transcript
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub transcript {
  return $_[0]->container->transcript;
}


=head2 get_other_Haplotypes

  Example    : my @ths = $th->get_other_Haplotypes()
  Description: Get the partner TranscriptHaplotypes to this one i.e. 
               get all CDSHaplotypes that translate to this ProteinHaplotype OR
               get the ProteinHaplotype representing the translation of this CDSHaplotype
  Returntype : arrayref of Bio::EnsEMBL::Variation::TranscriptHaplotype
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_other_Haplotypes {
  my $self = shift;
  return [map {$self->container->_get_TranscriptHaplotype_by_hex($_)} @{$self->_other_hexes}];
}


=head2 count

  Example    : my $count = $th->count()
  Description: Counts the number of times this haplotype has been observed
               within the individuals defined in the container
  Returntype : int
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub count {
  my $self = shift;
  
  if(!exists($self->{count})) {
    $self->{count} = $self->container->{_counts}->{$self->_hex};
  }
  
  return $self->{count};
}


=head2 frequency

  Example    : my $frequency = $th->frequency()
  Description: Get the observed frequency of this haplotype amongst
               those defined in the container
  Returntype : float
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub frequency {
  my $self = shift;
  
  if(!exists($self->{frequency})) {
    my $total = $self->container->total_haplotype_count;
    my $count = $self->count;
    $self->{frequency} = $count / $total;
  }
  
  return $self->{frequency};
}


=head2 get_all_population_counts

  Example    : my %counts = %{$th->get_all_population_counts()}
  Description: Counts the number of times this haplotype has been observed
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
    
    foreach my $sample(keys %{$self->{samples}}) {
      $counts->{$_} += $self->{samples}->{$sample} for keys %{$sample_pop_hash->{$sample}};
    }
    
    $self->{population_counts} = $counts;
  }
  
  return $self->{population_counts};
}


=head2 get_all_population_frequencies

  Example    : my %freqs = %{$th->get_all_population_frequencies()}
  Description: Gets the observed frequency of this haplotype in each population
               defined in the container
  Returntype : hashref e.g. { $population_name => $freq }
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_population_frequencies {
  my $self = shift;
  
  if(!exists($self->{population_frequencies})) {
    my $totals = $self->container->total_population_counts;
    my $counts = $self->get_all_population_counts;
    
    my %freqs = map {$_ => $counts->{$_} / $totals->{$_}} keys %$counts;
    
    $self->{population_frequencies} = \%freqs;
  }
  
  return $self->{population_frequencies};
}

sub get_aligned_sequences {
  my $self = shift;

  if(!exists($self->{aligned_sequences})) {

    my $aln = $self->_get_SimpleAlign_obj();

    # get sequences
    my @aligned = $aln->each_seq;

    # pad
    my $longest = (sort {$a <=> $b} map {$_->length} @aligned)[-1];
    my @padded = map {$_->seq.('-' x ($longest - $_->length))} @aligned;

    $self->{aligned_sequences} = \@padded;
  }

  return $self->{aligned_sequences};
}

sub get_formatted_alignment {
  my $self = shift;
  my $format = shift || 'clustalw';

  if(
    !exists($self->{formatted_alignment}) ||
    ($self->{alignment_format} && ($self->{alignment_format} ne '$format'))
  ) {
    my $aln_out;
    open(my $fh, ">", \$aln_out);
    Bio::AlignIO->new(-format => $format, -fh => $fh)->write_aln($self->_get_SimpleAlign_obj);

    $self->{formatted_alignment} = $aln_out;
    $self->{alignment_format} = 'format';
  }

  return $self->{formatted_alignment};
}

sub _get_SimpleAlign_obj {
  my $self = shift;

  if(!exists($self->{_SimpleAlign_obj})) {

    # if there are any indels we need to align
    if($self->has_indel) {

      # try and use bioperl-ext packages
      eval q{use Bio::Tools::dpAlign};

      # fall back to slow pure perl NW algorithm from Bio::Ensembl::Variation::Utils::Sequence
      if($@) {
        $self->{_SimpleAlign_obj} = $self->_create_SimpleAlign_from_sequence_pair(align_seqs($self->reference_seq, $self->seq));
      }

      # ok to use fast dpAlign
      else {

        # create Bio::Seq objects to feed aligner
        my $s1 = Bio::Seq->new(-id => "REF", -seq => $self->reference_seq);
        my $s2 = Bio::Seq->new(-id => "ALT", -seq => $self->seq);

        # get factory
        my $factory;

        eval {$factory = Bio::Tools::dpAlign->new(-alg => 2)};

        throw($@) if $@;

        # create alignment
        $self->{_SimpleAlign_obj} = $factory->pairwise_alignment($s1, $s2);
      }
    }

    else {
      $self->{_SimpleAlign_obj} = $self->_create_SimpleAlign_from_sequence_pair($self->reference_seq, $self->seq);
    }
  }

  return $self->{_SimpleAlign_obj};
}

sub _create_SimpleAlign_from_sequence_pair {
  my ($self, $s1, $s2) = @_;

  my $aln = Bio::SimpleAlign->new();
  $aln->add_seq(Bio::LocatableSeq->new(-id => "REF", -seq => $s1));
  $aln->add_seq(Bio::LocatableSeq->new(-id => "ALT", -seq => $s2));

  return $aln;
}

## Get/set the hex string uniquely identifying this haplotype
sub _hex {
  my $self = shift;
  $self->{hex} = shift if @_;
  return $self->{hex};
}

## Get/set an arrayref of hexes identifying the partner haplotypes to this one
## i.e. all possible CDSHaplotypes for a ProteinHaplotype
## or the ProteinHaplotype for a CDSHaplotype
sub _other_hexes {
  my $self = shift;
  $self->{other_hexes} = shift if @_;
  return [keys %{$self->{other_hexes}}];
}

## Get differences between this TranscriptHaplotype and the reference sequence
## Uses align_seqs (slow!) to do alignment if there's an indel
sub _get_raw_diffs {
  my $self = shift;
  
  if(!exists($self->{_raw_diffs})) {

    my ($al1, $al2) = @{$self->get_aligned_sequences};
    my @diffs;
    
    if($al1 ne $al2) {

      # if there was a stop introduced, then $al2 will be shorter
      # pad it out with '-' to make seqs same length
      $al2 .= '-' x (length($al1) - length($al2));
      
      ## this code finds diffs between the sequences
      ## copied verbatim from http://www.perlmonks.org/?node_id=882593
      ## can probably be XS'd in future, see http://www.perlmonks.org/?node_id=882590
      my ($m, @pos, $c1, @p1, $c2, @p2);
      $m  = $al1 ^ $al2;
      push @pos, pos( $m ) while $m =~ m{(?=([^\x00]))}g;
      $m  =~ tr{[\x01-\xfe]}{\xff};
      $c1 = $al1 & $m;
      @p1 = $c1 =~ m{([^\x00])}g;
      $c2 = $al2 & $m;
      @p2 = $c2 =~ m{([^\x00])}g;
      
      ## this bit then joins together consecutive mismatches
      ## but only if they are the same "type"
      ## i.e. join consecutive - or [ACGT] characters, but
      ## don't join - to A
      ## hopefully this doesn't make it too slow after the efforts above
      my ($p, $pp, $a1, $a2);
      for my $i(0..$#pos) {
        $p = $pos[$i];
        
        # only check join if this isn't the first pos
        if(defined($pp)) {
          
          # check positions are consecutive
          # and consecutive alleles on each string are of same type
          if(
            ($pp == $p - 1) &&
            (($p1[$i] eq '-') == ($p1[$i-1] eq '-')) &&
            (($p2[$i] eq '-') == ($p2[$i-1] eq '-'))
          ) {
            # extend a1 and a2
            $a1 .= $p1[$i];
            $a2 .= $p2[$i];
          }
          
          # if not, create a new diff
          else {
            push @diffs, $self->_create_diff($a1, $a2, $pp);
            
            # re-initiate a1 and a2
            $a1 = $p1[$i];
            $a2 = $p2[$i];
          }
        }
        
        # initiate a1 and a2
        else {
          $a1 = $p1[$i];
          $a2 = $p2[$i];
        }
        
        # record previous position
        $pp = $p;
      }
      
      # there will be a diff left over
      push @diffs, $self->_create_diff($a1, $a2, $pp);
    }
    
    $self->{_raw_diffs} = \@diffs;
  }
  
  return $self->{_raw_diffs};
}

## creates a diff string from a pair of "alleles" and a position
sub _create_diff {
  my ($self, $a1, $a2, $pp) = @_;
  
  my $pos = ($pp - length($a1)) + 2;
  
  # insertion
  if($a1 =~ /\-/) {
    $a2 = '{'.length($a2).'}' if length($a2) > 3;

    return $pos."ins".$a2;
  }

  # deletion
  elsif($a2 =~ /\-/) {
    $a1 = '{'.length($a1).'}' if length($a1) > 3;

    return $pos."del".$a1;
  }

  # substitution
  else {
    return $pos."$a1\>$a2";
  }
}

## Convert this object to a hash that can be written as JSON.
## Basically just deletes "private" keys starting with "_"
sub TO_JSON {
  my $self = shift;
  
  # make a hash copy of self
  my %copy = %{$self};
  
  # delete keys starting with _
  delete $copy{$_} for grep {$_ =~ /^\_/} keys %copy;
  
  return \%copy;
}

1;
