=head1 LICENSE

Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

A helper class for representing a transcript sequence modified by individual
genotypes. Not to be used directly.

=cut

package Bio::EnsEMBL::Variation::TranscriptHaplotype;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::Utils::Sequence qw(align_seqs);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  my ($hex, $type, $seq, $indel, $container) = rearrange([qw(HEX TYPE SEQ INDEL CONTAINER)], @_);
  
  my $self = {
    type => $type,
    seq => $seq,
    has_indel => $indel,
    hex => $hex,
    _container => $container
  };
  
  bless($self, $class);
  
  return $self;
}

sub container {
  my $self = shift;
  $self->{_container} = shift if @_;
  return $self->{_container};
}

sub type {
  my $self = shift;
  $self->{type} = shift if @_;
  return $self->{type};
}

sub seq {
  my $self = shift;
  $self->{seq} = shift if @_;
  return $self->{seq};
}

sub has_indel {
  my $self = shift;
  $self->{has_indel} = shift if @_;
  return $self->{has_indel};
}

sub hex {
  my $self = shift;
  $self->{hex} = shift if @_;
  return $self->{hex};
}

sub other_hexes {
  my $self = shift;
  $self->{other_hexes} = shift if @_;
  return $self->{other_hexes};
}

sub transcript {
  return $_[0]->container->transcript;
}

sub get_other_Haplotypes {
  my $self = shift;
  return [map {$self->container->get_TranscriptHaplotype_by_hex($_)} @{$self->other_hexes}];
}

sub count {
  my $self = shift;
  
  if(!exists($self->{count})) {
    $self->{count} = $self->container->{_counts}->{$self->hex};
  }
  
  return $self->{count};
}

sub frequency {
  my $self = shift;
  
  if(!exists($self->{frequency})) {
    my $total = $self->container->total_haplotype_count;
    my $count = $self->count;
    $self->{frequency} = $count / $total;
  }
  
  return $self->{frequency};
}

sub get_all_population_counts {
  my $self = shift;
  
  if(!exists($self->{population_counts})) {
    my $ind_pop_hash = $self->container->_get_individual_population_hash();
    my $counts = {};
    
    foreach my $ind(keys %{$self->{individuals}}) {
      $counts->{$_} += $self->{individuals}->{$ind} for keys %{$ind_pop_hash->{$ind}};
    }
    
    $self->{population_counts} = $counts;
  }
  
  return $self->{population_counts};
}

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

sub name {
  my $self = shift;
  
  if(!exists($self->{name})) {
    $self->{name} = sprintf("%s:%s", $self->transcript->stable_id, join(",", @{$self->_get_raw_diffs}) || 'REF');
  }
  
  return $self->{name};
}

## Get differences between this TranscriptHaplotype and the reference sequence
## Uses align_seqs (slow!) to do alignment if there's an indel
sub _get_raw_diffs {
  my $self = shift;
  
  if(!exists($self->{_raw_diffs})) {
    my $s1 = $self->transcript->{$self->type};
    my $s2 = $self->seq;
    my $indel = $self->has_indel;
    
    my @diffs;
    
    if($s1 ne $s2) {
      my ($al1, $al2) = $indel ? @{align_seqs($s1, $s2)} : ($s1, $s2);
      
      # if there was a stop introduced, then $al2 will be shorter
      # pad it out with '-' to make seqs same length
      $al2 .= '-' x (length($al1) - length($al2));  
      
      for(my $i=0; $i<length($al1); $i++) {
        my $len = 1;
        
        my $a1 = substr($al1, $i, $len);
        my $a2 = substr($al2, $i, $len);
        
        # try to extend if indel
        while(($a1 =~ /^\-+$/ || $a2 =~ /^\-+$/) && $i + $len <= length($al1)) {
          $len++;
          
          $a1 = substr($al1, $i, $len);
          $a2 = substr($al2, $i, $len);
        }
        
        # wind it back one if we extended as we will have gone one pos further
        if($len > 1) {
          $len--;
          $a1 = substr($al1, $i, $len);
          $a2 = substr($al2, $i, $len);
        }
        
        if($a1 ne $a2) {
          
          # insertion
          if($a1 =~ /\-/) {
            $a2 = '{'.length($a2).'}' if length($a2) > 3;
            
            push @diffs, ($i + 1)."ins".$a2;
          }
          
          # deletion
          elsif($a2 =~ /\-/) {
            $a1 = '{'.length($a1).'}' if length($a1) > 3;
            
            push @diffs, ($i + 1)."del".$a1;
          }
          
          # substitution
          else {
            push @diffs, ($i + 1)."$a1\>$a2";
          }
        }
        
        $i += $len - 1;
      }
    }
    
    $self->{_raw_diffs} = \@diffs;
  }
  
  return $self->{_raw_diffs};
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