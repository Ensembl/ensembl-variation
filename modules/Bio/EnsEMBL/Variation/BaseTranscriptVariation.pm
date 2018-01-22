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

Bio::EnsEMBL::Variation::BaseTranscriptVariation

=head1 SYNOPSIS

    use Bio::EnsEMBL::Variation::BaseTranscriptVariation;

=head1 DESCRIPTION

A helper class for representing an overlap of a Transcript and a
Variation (either sequence or structural). Should not be invoked directly.

=cut

package Bio::EnsEMBL::Variation::BaseTranscriptVariation;

use strict;
use warnings;

use Digest::MD5 qw(md5_hex);

use Bio::EnsEMBL::Utils::Scalar qw(assert_ref check_ref);
use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(overlap within_cds);

use base qw(Bio::EnsEMBL::Variation::VariationFeatureOverlap);

our $CAN_USE_INTERVAL_TREE;

BEGIN {
  if (eval { require Set::IntervalTree; 1 }) {
    $CAN_USE_INTERVAL_TREE = 1;
  }
  else {
    $CAN_USE_INTERVAL_TREE = 0;
  }
}


=head2 transcript_stable_id

  Description: Returns the stable_id of the associated Transcript
  Returntype : string
  Exceptions : none
  Status     : Stable

=cut

sub transcript_stable_id {
    my $self = shift;
    return $self->SUPER::_feature_stable_id(@_);
}

=head2 transcript

  Arg [1]    : (optional) Bio::EnsEMBL::Transcript
  Description: Get/set the associated Transcript
  Returntype : Bio::EnsEMBL::Transcript
  Exceptions : throws if argument is wrong type
  Status     : Stable

=cut

sub transcript {
  my ($self, $transcript) = @_;

  if($transcript) {
    assert_ref($transcript, 'Bio::EnsEMBL::Transcript') if $Bio::EnsEMBL::Utils::Scalar::ASSERTIONS && $transcript;
    delete $self->{_cached_transcript};
  }

  $self->{_cached_transcript} ||= $self->SUPER::feature($transcript, 'Transcript');

  return $self->{_cached_transcript};
}

=head2 feature

  Arg [1]    : (optional) Bio::EnsEMBL::Transcript
  Description: Get/set the associated Transcript (overriding the superclass feature method)
  Returntype : Bio::EnsEMBL::Transcript
  Exceptions : throws if argument is wrong type
  Status     : Stable

=cut

sub feature {
    my $self = shift;
    return $self->transcript(@_);
}

=head2 cdna_start

  Arg [1]    : (optional) int $start
  Example    : $cdna_start = $tv->cdna_start;
  Description: Getter/Setter for the start position of this variation on the
               transcript in cDNA coordinates.
  Returntype : int
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut

sub cdna_start {
    my ($self, $cdna_start) = @_;
    
    $self->{cdna_start} = $cdna_start if defined $cdna_start;
    
    unless (exists $self->{cdna_start}) {
        my $cdna_coords = $self->cdna_coords;
        
        my ($first, $last) = ($cdna_coords->[0], $cdna_coords->[-1]);
        
        $self->{cdna_start} = $first->isa('Bio::EnsEMBL::Mapper::Gap') ? undef : $first->start;
        $self->{cdna_end}   = $last->isa('Bio::EnsEMBL::Mapper::Gap') ? undef : $last->end;
    }
    
    return $self->{cdna_start};
}

=head2 cdna_end

  Arg [1]    : (optional) int $end
  Example    : $cdna_end = $tv->cdna_end;
  Description: Getter/Setter for the end position of this variation on the
               transcript in cDNA coordinates.
  Returntype : int
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut

sub cdna_end {
    my ($self, $cdna_end) = @_;
    
    $self->{cdna_end} = $cdna_end if defined $cdna_end;
    
    # call cdna_start to calculate the start and end
    $self->cdna_start unless exists $self->{cdna_end};
    
    return $self->{cdna_end};
}

=head2 cds_start

  Arg [1]    : (optional) int $start
  Example    : $cds_start = $tv->cds_start;
  Description: Getter/Setter for the start position of this variation on the
               transcript in CDS coordinates.
  Returntype : int
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut

sub cds_start {
    my ($self, $cds_start) = @_;
    
    $self->{cds_start} = $cds_start if defined $cds_start;
    
    unless (exists $self->{cds_start}) {
        my $cds_coords = $self->cds_coords;
        
        my ($first, $last) = ($cds_coords->[0], $cds_coords->[-1]);
        my $exon_phase = $self->transcript->start_Exon->phase;
        
        $self->{cds_start} = $first->isa('Bio::EnsEMBL::Mapper::Gap') ? undef : $first->start + ($exon_phase > 0 ? $exon_phase : 0);
        $self->{cds_end}   = $last->isa('Bio::EnsEMBL::Mapper::Gap') ? undef : $last->end + ($exon_phase > 0 ? $exon_phase : 0);
    }
    
    return $self->{cds_start};
}

=head2 cds_end

  Arg [1]    : (optional) int $end
  Example    : $cds_end = $tv->cds_end;
  Description: Getter/Setter for the end position of this variation on the
               transcript in CDS coordinates.
  Returntype : int
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut

sub cds_end {
    my ($self, $cds_end) = @_;
    
    $self->{cds_end} = $cds_end if defined $cds_end;
    
    # call cds_start to calculate the start and end
    $self->cds_start unless exists $self->{cds_end};
    
    return $self->{cds_end};
}

=head2 translation_start

  Arg [1]    : (optional) int $start
  Example    : $translation_start = $tv->translation_start;
  Description: Getter/Setter for the start position of this variation on the
               transcript in peptide coordinates.
  Returntype : int
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut

sub translation_start {
    my ($self, $translation_start) = @_;
    
    $self->{translation_start} = $translation_start if defined $translation_start;
    
    unless (exists $self->{translation_start}) {
        my $translation_coords = $self->translation_coords;
        
        my ($first, $last) = ($translation_coords->[0], $translation_coords->[-1]);
        
        $self->{translation_start} = $first->isa('Bio::EnsEMBL::Mapper::Gap') ? undef : $first->start;
        $self->{translation_end}   = $last->isa('Bio::EnsEMBL::Mapper::Gap') ? undef : $last->end;
    }

    return $self->{translation_start};
}


=head2 translation_end

  Arg [1]    : (optional) int $end
  Example    : $transaltion_end = $tv->translation_end;
  Description: Getter/Setter for the end position of this variation on the
               transcript in peptide coordinates.
  Returntype : int
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut

sub translation_end {
    my ($self, $translation_end) = @_;
    
    $self->{translation_end} = $translation_end if defined $translation_end;
    
    # call translation_start to calculate the start and end
    $self->translation_start unless exists $self->{translation_end};
    
    return $self->{translation_end};
}

=head2 cdna_coords

  Description: Use the TranscriptMapper to calculate the cDNA 
               coordinates of this variation
  Returntype : listref of Bio::EnsEMBL::Coordinate and Bio::EnsEMBL::Gap objects
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut

sub cdna_coords {
    my $self = shift;
    
    unless ($self->{_cdna_coords}) {
        my $vf   = $self->base_variation_feature;
        my $tran = $self->transcript;
        $self->{_cdna_coords} = [ $self->_mapper->genomic2cdna($vf->seq_region_start, $vf->seq_region_end, $tran->strand) ];
    }
    
    return $self->{_cdna_coords};
}

=head2 cds_coords

  Description: Use the TranscriptMapper to calculate the CDS 
               coordinates of this variation
  Returntype : listref of Bio::EnsEMBL::Coordinate and Bio::EnsEMBL::Gap objects
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut

sub cds_coords {
    my $self = shift;
    
    unless ($self->{_cds_coords}) {
        my $vf   = $self->base_variation_feature;
        my $tran = $self->transcript;
        $self->{_cds_coords} = [ $self->_mapper->genomic2cds($vf->seq_region_start, $vf->seq_region_end, $tran->strand) ];
    }
    
    return $self->{_cds_coords};
}

=head2 translation_coords

  Description: Use the TranscriptMapper to calculate the peptide
               coordinates of this variation
  Returntype : listref of Bio::EnsEMBL::Coordinate and Bio::EnsEMBL::Gap objects
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut

sub translation_coords {
    my $self = shift;
    
    unless ($self->{_translation_coords}) {
        my $vf   = $self->base_variation_feature;
        my $tran = $self->transcript; 
        $self->{_translation_coords} = [ $self->_mapper->genomic2pep($vf->seq_region_start, $vf->seq_region_end, $tran->strand) ];
    }
    
    return $self->{_translation_coords};
}


=head2 distance_to_transcript

  Arg [1]    : (optional) int $distance
  Example    : $distance = $tv->distance_to_transcript;
  Description: Getter/Setter for the distance of this variant to the transcript.
               This is the shortest distance between variant start/end and
               transcript start/end, so if a variant falls 5' of a transcript
               on the forward strand this distance will be that between the
               variant end and the transcript start; if it falls 3' it will be
               the distance between the variant start and the transcript end.
  Returntype : int
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut

sub distance_to_transcript {
    my ($self, $distance) = @_;
    
    $self->{distance_to_transcript} = $distance if defined $distance;
    
    unless (exists $self->{distance_to_transcript}) {
        my $vf = $self->base_variation_feature;
        my $tr = $self->transcript;
        
        my @dists = (
            $vf->{start} - $tr->{start},
            $vf->{start} - $tr->{end},
            $vf->{end} - $tr->{start},
            $vf->{end} - $tr->{end}
        );
        
        # make positive if <0 and sort
        @dists = sort {$a <=> $b} map {abs($_)} @dists;
        
        $self->{distance_to_transcript} = $dists[0];
    }
    
    return $self->{distance_to_transcript};
}

=head2 get_overlapping_ProteinFeatures

  Description: Find any ProteinFeatures (e.g. pfam or interpro domains etc.) that
               the associated variation feature lies in
  Returntype : listref of Bio::EnsEMBL::ProteinFeatures (possibly empty)
  Exceptions : None
  Caller     : general
  Status     : At Risk

=cut

sub get_overlapping_ProteinFeatures {
    my $self = shift;

    unless (exists $self->{_protein_features}) {

        $self->{_protein_features } = [];

        my $tr = $self->transcript;
        my $tl = $tr->translation;

        if (defined $tl) {
            
            my $tl_start = $self->translation_start;
            my $tl_end   = $self->translation_end;

            if (defined $tl_start && defined $tl_end) {
                unless(exists $tr->{_variation_effect_feature_cache}->{protein_features}) {
                    my @feats = @{$tl->get_all_ProteinFeatures || []};
                    $tr->{_variation_effect_feature_cache}->{protein_features} = \@feats;
                }
              
                for my $feat (@{ $tr->{_variation_effect_feature_cache}->{protein_features} }) {
                    if (overlap($feat->{start}, $feat->{end}, $tl_start, $tl_end)) { 
                        push @{ $self->{_protein_features} }, $feat;
                    }
                }
            }
        }
    }

    return $self->{_protein_features};
}

=head2 affects_cds

  Description: Check if any of this TranscriptVariation's alleles lie within
               the CDS of the Transcript
  Returntype : boolean
  Exceptions : None
  Caller     : general
  Status     : At Risk

=cut

sub affects_cds {
    my $self = shift;
    return scalar grep { within_cds($_) } @{ $self->get_all_alternate_BaseVariationFeatureOverlapAlleles }; 
}

=head2 exon_number

  Description: Identify which exon(s) this variant falls in   
  Returntype : '/'-separated string containing the exon number(s) and the total 
               number of exons in this transcript, or undef if this variant 
               does not fall in any exons
  Exceptions : None
  Caller     : general
  Status     : At Risk

=cut

sub exon_number {
  my $self = shift;

  if(!exists($self->{_exon_number})) {
    my $total = scalar @{$self->_exons};

    my $tr = $self->transcript;
    my $vf = $self->base_variation_feature;
    my $vf_start = $vf->{start};
    my $vf_end   = $vf->{end};

    my @numbers =
      map {$tr->{_variation_effect_feature_cache}->{_exon_numbers}->{sprintf('%s', $_)}}
      grep {overlap($vf_start, $vf_end, $_->{start}, $_->{end})}
      @{$self->_overlapped_exons};

    my $number_string = undef;

    if(@numbers) {
      if(scalar @numbers > 1) {
        @numbers = sort {$a <=> $b} @numbers;
        $number_string = $numbers[0].'-'.$numbers[-1];
      }
      else {
        $number_string = $numbers[0];
      }

      $number_string .= '/'.$total;
    }

    $self->{_exon_number} = $number_string;
  }

  return $self->{_exon_number};
}

=head2 intron_number
  
  Description: Identify which intron(s) this variant falls in   
  Returntype : '/'-separated string containing the intron number(s) and the total 
               number of introns in this transcript, or undef if this variant 
               does not fall in any introns
  Exceptions : None
  Caller     : general
  Status     : At Risk

=cut

sub intron_number {
  my $self = shift;

  if(!exists($self->{_intron_number})) {
    my $total = scalar @{$self->_introns};

    my $tr = $self->transcript;
    my $vf = $self->base_variation_feature;
    my $vf_start = $vf->{start};
    my $vf_end   = $vf->{end};

    my @numbers =
      map {$tr->{_variation_effect_feature_cache}->{_intron_numbers}->{sprintf('%s', $_)}}
      grep {overlap($vf_start, $vf_end, $_->{start}, $_->{end})}
      @{$self->_overlapped_introns};

    my $number_string = undef;

    if(@numbers) {
      if(scalar @numbers > 1) {
        @numbers = sort {$a <=> $b} @numbers;
        $number_string = $numbers[0].'-'.$numbers[-1];
      }
      else {
        $number_string = $numbers[0];
      }

      $number_string .= '/'.$total;
    }

    $self->{_intron_number} = $number_string;
  }

  return $self->{_intron_number};
}

# NB: the methods below all cache their data in the associated transcript itself, this
# gives a significant speed up when you are calculating the effect of all variations
# on a transcript, and means that the cache will be freed when the transcript itself
# is garbage collected rather than us having to maintain a transcript feature cache 
# ourselves

# gets all introns for the transcript, cached on transcript
sub _introns {
  my $self = shift;
  
  my $tran = $self->transcript;
  
  my $introns = $tran->{_variation_effect_feature_cache}->{introns} ||= $tran->get_all_Introns;
  
  return $introns;
}

# gets all exons for the transcript, cached on transcript
sub _exons {
  my $self = shift;

  my $tran = $self->transcript;

  my $exons = $tran->{_variation_effect_feature_cache}->{exons} ||= $tran->get_all_Exons;

  return $exons;
}

sub _sorted_exons {
  my $self = shift;

  my $tran = $self->transcript;

  my $exons = $tran->{_variation_effect_feature_cache}->{sorted_exons} ||= [sort {$a->start <=> $b->start} @{$self->_exons}];

  return $exons;
}

# gets all introns that overlap this VF, cached on $self
sub _overlapped_introns {
  my ($self, $min_vf, $max_vf) = @_;

  if(!exists($self->{_overlapped_introns})) {

    if(!($min_vf || $max_vf)) {
      my $vf = $self->base_variation_feature;
      ($min_vf, $max_vf) = sort {$a <=> $b} ($vf->{start}, $vf->{end});
    }

    if($CAN_USE_INTERVAL_TREE) {
      $self->{_overlapped_introns} = $self->_intron_interval_tree->fetch($min_vf - 1, $max_vf);
    }

    else {
      $self->_overlapped_introns_and_boundary_no_tree($min_vf, $max_vf);
    }
  }

  return $self->{_overlapped_introns};
}

# gets all introns whose boundary region overlaps this VF, cached on $self
sub _overlapped_introns_boundary {
  my ($self, $min_vf, $max_vf) = @_;

  if(!exists($self->{_overlapped_introns_boundary})) {

    if(!($min_vf || $max_vf)) {
      my $vf = $self->base_variation_feature;
      ($min_vf, $max_vf) = sort {$a <=> $b} ($vf->{start}, $vf->{end});
    }

    if($CAN_USE_INTERVAL_TREE) {
      $self->{_overlapped_introns_boundary} = $self->_intron_boundary_interval_tree->fetch($min_vf - 1, $max_vf);
    }

    else {
      $self->_overlapped_introns_and_boundary_no_tree($min_vf, $max_vf);
    }
  }

  return $self->{_overlapped_introns_boundary};
}

# gets all exons that overlap this VF, cached on $self
sub _overlapped_exons {
  my ($self, $min_vf, $max_vf) = @_;

  if(!exists($self->{_overlapped_exons})) {

    if(!($min_vf || $max_vf)) {
      my $vf = $self->base_variation_feature;
      ($min_vf, $max_vf) = sort {$a <=> $b} ($vf->{start}, $vf->{end});
    }

    # apply a "stretch" to the exon coordinates if we have a frameshift intron
    # the introns can be up to 12 bases long and if a variant falls in one
    # we actually want to call it exonic
    my $stretch = $self->transcript->{_variation_effect_feature_cache}->{_has_frameshift_intron} ? 12 : 0;

    $self->{_overlapped_exons} = 
      $CAN_USE_INTERVAL_TREE ?
      $self->_exon_interval_tree->fetch($min_vf - ($stretch + 1), $max_vf + $stretch) :
      $self->_overlapped_exons_no_tree($min_vf, $max_vf, $stretch);
  }

  return $self->{_overlapped_exons};
}

# slower method for looking up overlaps
# used when user does not have Set::IntervalTree installed
sub _overlapped_introns_and_boundary_no_tree {
  my ($self, $min_vf, $max_vf) = @_;

  my $tr = $self->transcript;
  my $tr_strand = $tr->strand();
  my $vf_3_prime_end = $tr_strand > 0 ? $max_vf : $min_vf;
  
  my (@introns, @boundary);
  my $intron_number = 0;
  my $num_cache = $tr->{_variation_effect_feature_cache}->{_intron_numbers} ||= {};
  
  foreach my $intron(@{$self->_introns}) {

    # log numbers
    $num_cache->{sprintf('%s', $intron)} = ++$intron_number;

    my ($intron_start, $intron_end) = ($intron->{start}, $intron->{end});

    # sometimes RefSeq transcripts can have overlapping exons
    # this leads to "introns" coming out with start > end
    next unless $intron_start <= $intron_end;

    # check within intron
    if(overlap($min_vf, $max_vf, $intron_start, $intron_end)) {
      push @introns, $intron;
    }

    # check intron boundary
    # add flank to capture splice_region_variants
    if(
      overlap($min_vf, $max_vf, $intron_start - 3, $intron_start + 7) ||
      overlap($min_vf, $max_vf, $intron_end - 7, $intron_end + 3)
    ) {
      push @boundary, $intron;
    }

    # cache whether this transcript has a frameshift intron
    # we use it later to assess exonic status
    if(abs($intron_end - $intron_start) <= 12) {
      $tr->{_variation_effect_feature_cache}->{_has_frameshift_intron} = 1;
      $intron->{_frameshift} = 1;
    }
    
    # also leave if we've gone beyond the bounds of the VF
    # subsequent introns won't be in range
    last if
      ( $tr_strand > 1 && $vf_3_prime_end < ($intron->{start} - 3)) ||
      ( $tr_strand < 1 && $vf_3_prime_end > ($intron->{end} + 3));
  }

  $self->{_overlapped_introns} = \@introns;
  $self->{_overlapped_introns_boundary} = \@boundary;
}

# slower method for looking up overlaps
# used when user does not have Set::IntervalTree installed
sub _overlapped_exons_no_tree {
  my ($self, $min_vf, $max_vf, $stretch) = @_;

  my $tr = $self->transcript;
  my $tr_strand = $tr->strand();
  my $vf_3_prime_end = $tr_strand > 0 ? $max_vf : $min_vf;

  my @exons;
  my $exon_number = 0;
  my $num_cache = $tr->{_variation_effect_feature_cache}->{_exon_numbers} ||= {};

  foreach my $exon(@{$self->_exons}) {

    # log numbers
    $num_cache->{sprintf('%s', $exon)} = ++$exon_number;

    # also leave if we've gone beyond the bounds of the VF
    # subsequent introns won't be in range
    last if
      ( $tr_strand > 1 && $vf_3_prime_end < $exon->{start} - $stretch) ||
      ( $tr_strand < 1 && $vf_3_prime_end > $exon->{end} + $stretch);

    if(overlap($min_vf, $max_vf, $exon->{start} - $stretch, $exon->{end} + $stretch)) {
      push @exons, $exon;
    }
  }

  return \@exons;
}


## interval tree methods
########################

# create intron interval tree for this transcript, cached on transcript
sub _intron_interval_tree {
  my $self = shift;

  my $tr = $self->transcript;
  
  $self->_create_intron_trees() unless exists($tr->{_variation_effect_feature_cache}->{_intron_interval_tree});

  return $tr->{_variation_effect_feature_cache}->{_intron_interval_tree};
}

# create intron boundary interval tree for this transcript, cached on transcript
sub _intron_boundary_interval_tree {
  my $self = shift;

  my $tr = $self->transcript;
  
  $self->_create_intron_trees() unless exists($tr->{_variation_effect_feature_cache}->{_intron_boundary_interval_tree});

  return $tr->{_variation_effect_feature_cache}->{_intron_boundary_interval_tree};
}

# this creates the intron and intron boundary trees
# not cached as it should only be called by _intron_interval_tree or _intron_boundary_interval_tree
sub _create_intron_trees {
  my $self = shift;

  my $tr = $self->transcript;

  if(!$CAN_USE_INTERVAL_TREE) {
    $tr->{_variation_effect_feature_cache}->{_intron_interval_tree} = undef;
    $tr->{_variation_effect_feature_cache}->{_intron_boundary_interval_tree} = undef;
    return;
  }

  my $intron_tree   = Set::IntervalTree->new();
  my $boundary_tree = Set::IntervalTree->new();

  my $intron_number = 0;
  my $num_cache = $tr->{_variation_effect_feature_cache}->{_intron_numbers} ||= {};

  for my $intron(@{$self->_introns}) {

    # log numbers
    $num_cache->{sprintf('%s', $intron)} = ++$intron_number;

    my ($intron_start, $intron_end) = ($intron->{start}, $intron->{end});

    # sometimes RefSeq transcripts can have overlapping exons
    # this leads to "introns" coming out with start > end
    next unless $intron_start <= $intron_end;

    # this is the actual plot
    #
    # ..........IIIIIIIIIIIIIIIIIIII..........
    # .......rrrrrrrrrrrrrrrrrrrrrrrrrr.......
    # .......bbbbbbbbbbb....bbbbbbbbbbb.......
    #
    # I is the intron itself
    # r is the region we care about overlapping for intron
    # b is the boundary region for splice_region_variant calls
    # starts adjusted by -1 to account for 0-based coords in interval tree

    $intron_tree->insert($intron, $intron_start - 4, $intron_end + 3);
    $boundary_tree->insert($intron, $intron_start - 4, $intron_start + 7);
    $boundary_tree->insert($intron, $intron_end - 8, $intron_end + 3);

    # cache this as it affects whether we should call something as overlapping an exon
    if(abs($intron_end - $intron_start) <= 12) {
      $tr->{_variation_effect_feature_cache}->{_has_frameshift_intron} = 1;
      $intron->{_frameshift} = 1;
    }
  }

  $tr->{_variation_effect_feature_cache}->{_intron_interval_tree} = $intron_tree;
  $tr->{_variation_effect_feature_cache}->{_intron_boundary_interval_tree} = $boundary_tree;
}

# create exon interval tree for this transcript, cached on transcript
sub _exon_interval_tree {
  my $self = shift;

  my $tr = $self->transcript;

  if(!exists($tr->{_variation_effect_feature_cache}->{_exon_interval_tree})) {

    if(!$CAN_USE_INTERVAL_TREE) {
      $tr->{_variation_effect_feature_cache}->{_exon_interval_tree} = undef;
    }

    else {
      my $tree = Set::IntervalTree->new();
      my $exon_number = 0;
      my $num_cache = $tr->{_variation_effect_feature_cache}->{_exon_numbers} ||= {};

      foreach my $exon(@{$self->_exons}) {

        # add it to the tree
        $tree->insert($exon, $exon->{start} - 1, $exon->{end});

        # log numbers
        $num_cache->{sprintf('%s', $exon)} = ++$exon_number;
      }
      $tr->{_variation_effect_feature_cache}->{_exon_interval_tree} = $tree;
    }
  }

  return $tr->{_variation_effect_feature_cache}->{_exon_interval_tree};
}

sub _mapper {
    my $self = shift;
    
    my $tran = $self->transcript;
    
    my $mapper = $tran->{_variation_effect_feature_cache}->{mapper} ||= $tran->get_TranscriptMapper;
    
    return $mapper;
}
sub _translateable_seq {
    my $self = shift;
    
    my $tran = $self->transcript;
    
    my $tran_seq = $tran->{_variation_effect_feature_cache}->{translateable_seq} ||= $tran->translateable_seq;
    
    return $tran_seq;
}

sub _five_prime_utr {
  return $_[0]->_generic_utr('five');
}

sub _three_prime_utr {
  return $_[0]->_generic_utr('three');
}

sub _generic_utr {
    my ($self, $three_five) = @_;
    
    my $tran = $self->transcript;
    my $key = $three_five.'_prime_utr';

    unless(exists($tran->{_variation_effect_feature_cache}->{$key})) {
        
        # transfer to feature slice so we don't subseq whole chromosome
        my $transferred = $tran->transfer($tran->feature_Slice());
        $tran->{_variation_effect_feature_cache}->{$key} = $transferred->$key();
    }
    
    return $tran->{_variation_effect_feature_cache}->{$key};
}
  
sub _peptide {
    my $self = shift;
    
    my $tran = $self->transcript;
    
    my $peptide = $tran->{_variation_effect_feature_cache}->{peptide};
    
    unless ($peptide) {
        my $translation = $tran->translate;
        $peptide = $translation ? $translation->seq : undef;
        $tran->{_variation_effect_feature_cache}->{peptide} = $peptide;
    }
    
    return $peptide;
}

sub _translation_md5 {
    my $self = shift;

    my $tran = $self->transcript;
    
    unless (exists $tran->{_variation_effect_feature_cache}->{translation_md5}) {
        $tran->{_variation_effect_feature_cache}->{translation_md5} = 
            $self->_peptide ? md5_hex($self->_peptide) : undef;
    }

    return $tran->{_variation_effect_feature_cache}->{translation_md5};
}

sub _codon_table {
    my $self = shift;
    
    my $tran = $self->transcript;
    
    my $codon_table = $tran->{_variation_effect_feature_cache}->{codon_table};
    
    unless ($codon_table) {
        # for mithocondrial dna we need to to use a different codon table
        my $attrib = $tran->slice->get_all_Attributes('codon_table')->[0]; 
        
        # default to the vertebrate codon table which is denoted as 1
        $codon_table = $attrib ? $attrib->value : 1;
        
        $tran->{_variation_effect_feature_cache}->{codon_table} = $codon_table
    }
    
    return $codon_table;
}

sub _seq_edits {
  my $self = shift;
  
  my $tr = $self->transcript;
  
  unless(exists $tr->{_variation_effect_feature_cache}->{seq_edits}) {
    
    if($tr && ($self->adaptor->{db} || ($tr->adaptor && $tr->adaptor->{db}))) {
      my $tl = $tr->translation;
    
      if($tl) {
        $tr->{_variation_effect_feature_cache}->{seq_edits} = $tl->get_all_SeqEdits();
      }
    }
    
    $tr->{_variation_effect_feature_cache}->{seq_edits} ||= [];
  }
  
  return $tr->{_variation_effect_feature_cache}->{seq_edits};
}

1;