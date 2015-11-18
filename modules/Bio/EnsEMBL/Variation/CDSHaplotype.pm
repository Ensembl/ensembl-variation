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

Bio::EnsEMBL::Variation::CDSHaplotype

=head1 SYNOPSIS

  use Bio::EnsEMBL::Variation::CDSHaplotype;

=head1 DESCRIPTION

A class for representing a transcript's CDS sequence modified by sample
genotypes.

=cut

package Bio::EnsEMBL::Variation::CDSHaplotype;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::Utils::Sequence qw(align_seqs);

use Bio::EnsEMBL::Variation::TranscriptHaplotype;

use base qw(Bio::EnsEMBL::Variation::TranscriptHaplotype);

=head2 new

  Arg [-CONTAINER]:  Bio::EnsEMBL::Variation::TranscriptHaplotypeContainer
  Arg [-SEQ]:        string
  Arg [-HEX]:        string
  Arg [-INDEL]:      bool

  Example    : my $ch = Bio::EnsEMBL::Variation::CDSHaplotype->new(
                  -CONTAINER => $container,
                  -SEQ       => $seq,
                  -HEX       => $hex,
                  -INDEL     => $indel
               );

  Description: Constructor.  Instantiates a new CDSHaplotype object.
  Returntype : Bio::EnsEMBL::Variation::DBSQL::CDSHaplotype
  Exceptions : none
  Caller     : internal
  Status     : Stable

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  my %args = @_;
  
  my $self = $class->SUPER::new(%args);
  bless($self, $class);

  $self->type('cds');
  
  return $self;
}


=head2 get_ProteinHaplotype

  Example    : my $ph = $ch->get_ProteinHaplotype()
  Description: Get the ProteinHaplotype representing the translation of this CDSHaplotype
  Returntype : Bio::EnsEMBL::Variation::ProteinHaplotype
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_ProteinHaplotype {
  return $_[0]->get_other_Haplotypes->[0];
}


=head2 reference_seq

  Example    : my $ref_seq = $ph->reference_seq
  Description: Get the reference CDS sequence
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub reference_seq {
  return $_[0]->transcript->{cds};
}


=head2 get_all_diffs

  Example    : my @diffs = @{$ch->get_all_diffs}
  Description: Get a list of differences to the reference. Each difference is a
               hashref containing a string 'diff' representing a change
               e.g. 25A>T represents a change of "A" to "T" at position 25. The
               hashref may also contain a pointer to the VariationFeature object
               that contributes this change (if a single VariatioFeature can be
               identified)
  Returntype : arrayref of hashrefs
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_diffs {
  my $self = shift;
  
  if(!defined($self->{diffs})) {
    my @diffs = ();
    
    my $raw_to_vf = $self->container->_get_raw_diff_to_VariationFeature_hash();
    
    foreach my $raw(@{$self->_get_raw_diffs()}) {
      my $diff = {
        'diff' => $raw, 
      };
      
      $diff->{variation_feature} = $raw_to_vf->{$raw} if $raw_to_vf->{$raw};
      
      push @diffs, $diff;
    }
    
    $self->{diffs} = \@diffs;
  }
  
  return $self->{diffs};
}

sub _reference_name {
  return $_[0]->transcript->stable_id;
}

sub TO_JSON {
  my $self = shift;
  
  # make a hash copy of self
  my %copy = %{$self};
  
  # delete keys starting with _
  delete $copy{$_} for grep {$_ =~ /^\_/} keys %copy;
  
  # change VF links to VF name
  if($self->{diffs}) {
    $copy{diffs} = [];

    foreach my $diff(@{$self->{diffs} || []}) {
      my %diff_copy = %{$diff};

      if($diff->{variation_feature}) {
        $diff_copy{variation_feature} = $diff->{variation_feature}->variation_name;
        $diff_copy{variation_feature_id} = $diff->{variation_feature}->dbID;
      }
      
      push @{$copy{diffs}}, \%diff_copy;
    }
  }
  
  return \%copy;
}

1;
