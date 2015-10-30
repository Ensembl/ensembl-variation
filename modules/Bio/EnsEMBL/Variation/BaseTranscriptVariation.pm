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
use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(overlap within_cds _intron_overlap);

use base qw(Bio::EnsEMBL::Variation::VariationFeatureOverlap);

our $CAN_USE_INTERVAL_TREE;

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
    assert_ref($transcript, 'Bio::EnsEMBL::Transcript') if $transcript;
    delete $self->{_cached_transcript};
  }

  $self->{_cached_transcript} ||= $self->SUPER::feature($transcript, 'Transcript');
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
        $self->{_cdna_coords} = [ $self->_mapper->genomic2cdna($vf->start, $vf->end, $tran->strand) ];
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
        $self->{_cds_coords} = [ $self->_mapper->genomic2cds($vf->start, $vf->end, $tran->strand) ];
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
        $self->{_translation_coords} = [ $self->_mapper->genomic2pep($vf->start, $vf->end, $tran->strand) ];
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
            $vf->start - $tr->start,
            $vf->start - $tr->end,
            $vf->end - $tr->start,
            $vf->end - $tr->end
        );
        
        # make positive if <0 and sort
        @dists = sort {$a <=> $b} map {$_ < 0 ? 0 - $_ : $_} @dists;
        
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
                    if (overlap($feat->start, $feat->end, $tl_start, $tl_end)) { 
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
    $self->_exon_intron_number unless exists $self->{exon_number};
    return $self->{exon_number};
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
    $self->_exon_intron_number unless exists $self->{intron_number};
    return $self->{intron_number};
}

sub _exon_intron_number {
    my $self = shift;

    # work out which exon or intron this variant falls in

    # ensure the keys exist so even if we don't fall in an exon 
    # or intron we'll only call this method once

    $self->{exon_number} = $self->{intron_number} = undef;

    my $vf = $self->base_variation_feature;    
    
    my $vf_start = $vf->start;
    my $vf_end   = $vf->end;

    my $strand = $self->transcript->strand;

    my $exons = $self->_exons;

    my $tot_exons = scalar(@$exons);

    my $exon_count = 0;

    my $prev_exon;
    
    my (@overlapped_exons, @overlapped_introns);

    for my $exon (@$exons) {

        $exon_count++;
        
        if (overlap($vf_start, $vf_end, $exon->start, $exon->end)) {
            push @overlapped_exons, $exon_count;
            #$self->{exon_number} = defined($self->{exon_number}) ? $self->{exon_number}.",".$exon_count : $exon_count;
        }

        if ($prev_exon) {
            my $intron_start = $strand == 1 ? $prev_exon->end + 1 : $exon->end + 1;
            my $intron_end   = $strand == 1 ? $exon->start - 1 : $prev_exon->start - 1;

            if ($prev_exon && overlap($vf_start, $vf_end, $intron_start, $intron_end)) {
                push @overlapped_introns, $exon_count - 1;
                #$self->{intron_number} = defined($self->{intron_number}) ? $self->{intron_number}.",".($exon_count - 1) : ($exon_count - 1);
            }
        }

        $prev_exon = $exon;
    }
    
    if(@overlapped_exons) {
        $self->{exon_number} = (scalar @overlapped_exons > 1 ? $overlapped_exons[0].'-'.$overlapped_exons[-1] : $overlapped_exons[0]).'/'.$tot_exons;
    }
    if(@overlapped_introns) {
        $self->{intron_number} = (scalar @overlapped_introns > 1 ? $overlapped_introns[0].'-'.$overlapped_introns[-1] : $overlapped_introns[0]).'/'.($tot_exons - 1);
    }
}

sub _intron_effects {
    my $self = shift;

    # internal method used by Bio::EnsEMBL::Variation::Utils::VariationEffect
    # when calculating various consequence types
    
    # this method is a major bottle neck in the effect calculation code so 
    # we cache results and use local variables instead of method calls where
    # possible to speed things up - caveat bug-fixer!
    
    unless ($self->{_intron_effects}) {
        
        my $vf = $self->base_variation_feature;
        
        my $intron_effects = {};
        
        my $found_effect = 0;
        
        my $vf_start = $vf->start;
        my $vf_end   = $vf->end;

        my $insertion = $vf_start == $vf_end+1;
        
        my $tr = $self->transcript();
        my $tr_strand = $tr->strand();
        
        my ($min, $max) = sort {$a <=> $b} ($vf_start, $vf_end);
        my $vf_3_prime_end = $tr_strand > 0 ? $max : $min;

        my @introns;

        my $tree = $self->_intron_interval_tree();
        if($tree) {
          @introns = @{$tree->fetch($min - 1, $max)};
        }
        else {
          @introns = @{$self->_introns};
        }

        for my $intron (@introns) {

            my $intron_start = $intron->start;
            my $intron_end   = $intron->end;
            
            # don't need these coord checks if we used the interval tree
            unless($tree) {
             
              # skip remainder if vf out of range
              last if
                ( $tr_strand > 1 && $vf_3_prime_end < $intron_start - 8) ||
                ( $tr_strand < 1 && $vf_3_prime_end > $intron_end + 8);
                
              # skip this one if vf out of range
              next unless overlap($min, $max, $intron_start - 8, $intron_end + 8);
            }
            
            # under various circumstances the genebuild process can introduce
            # artificial short (<= 12 nucleotide) introns into transcripts
            # (e.g. to deal with errors in the reference sequence etc.), we
            # don't want to categorise variations that fall in these introns
            # as intronic, or as any kind of splice variant
            
            my $frameshift_intron = ( abs($intron_end - $intron_start) <= 12 );
            
            if ($frameshift_intron) {
                if (overlap($vf_start, $vf_end, $intron_start, $intron_end)) {
                    $intron_effects->{within_frameshift_intron} = 1;
                    next;
                }
            }

            if (overlap($vf_start, $vf_end, $intron_start, $intron_start+1)) {
                $intron_effects->{start_splice_site} = 1;
            }
            
            if (overlap($vf_start, $vf_end, $intron_end-1, $intron_end)) {
                $intron_effects->{end_splice_site} = 1;
            }
            
            # we need to special case insertions between the donor and acceptor sites

            if (overlap($vf_start, $vf_end, $intron_start+2, $intron_end-2) or 
                ($insertion && ($vf_start == $intron_start+2 || $vf_end == $intron_end-2)) ) {
                $intron_effects->{intronic} = 1;
            }
            
            # the definition of splice_region (SO:0001630) is "within 1-3 bases 
            # of the exon or 3-8 bases of the intron", the intron start is the 
            # first base of the intron so we only need to add or subtract 7 from 
            # it to get the correct coordinate. We also need to special case 
            # insertions between the edge of an exon and a donor or acceptor site
            # and between a donor or acceptor site and the intron
            
            $intron_effects->{splice_region} = _intron_overlap($vf_start, $vf_end, $intron_start, $intron_end, $insertion);
        }
        
        $self->{_intron_effects} = $intron_effects;       
    }

    return $self->{_intron_effects};
}

# NB: the methods below all cache their data in the associated transcript itself, this
# gives a significant speed up when you are calculating the effect of all variations
# on a transcript, and means that the cache will be freed when the transcript itself
# is garbage collected rather than us having to maintain a transcript feature cache 
# ourselves

sub _introns {
    my $self = shift;
    
    my $tran = $self->transcript;
    
    my $introns = $tran->{_variation_effect_feature_cache}->{introns} ||= $tran->get_all_Introns;
    
    return $introns;
}

sub _exons {
    my $self = shift;

    my $tran = $self->transcript;

    my $exons = $tran->{_variation_effect_feature_cache}->{exons} ||= $tran->get_all_Exons;

    return $exons;
}

sub _can_use_interval_tree {
  my $self = shift;

  if(!defined($CAN_USE_INTERVAL_TREE)) {
    eval q{ use Set::IntervalTree; };
    $CAN_USE_INTERVAL_TREE = $@ ? 0 : 1;
  }

  return $CAN_USE_INTERVAL_TREE;
}

sub _intron_interval_tree {
  my $self = shift;

  my $tr = $self->transcript;
  
  $self->_create_intron_trees() unless exists($tr->{_variation_effect_feature_cache}->{_intron_interval_tree});

  return $tr->{_variation_effect_feature_cache}->{_intron_interval_tree};
}

sub _intron_boundary_interval_tree {
  my $self = shift;

  my $tr = $self->transcript;
  
  $self->_create_intron_trees() unless exists($tr->{_variation_effect_feature_cache}->{_intron_boundary_interval_tree});

  return $tr->{_variation_effect_feature_cache}->{_intron_boundary_interval_tree};
}

# this creates the intron and intron boundary trees
sub _create_intron_trees {
  my $self = shift;

  my $tr = $self->transcript;

  if(!$self->_can_use_interval_tree) {
    $tr->{_variation_effect_feature_cache}->{_intron_interval_tree} = undef;
    $tr->{_variation_effect_feature_cache}->{_intron_boundary_interval_tree} = undef;
    return;
  }

  my $intron_tree   = Set::IntervalTree->new();
  my $boundary_tree = Set::IntervalTree->new();

  for(@{$self->_introns}) {
    my ($intron_start, $intron_end) = ($_->start, $_->end);

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

    $intron_tree->insert($_, $intron_start - 4, $intron_end + 3);
    $boundary_tree->insert($_, $intron_start - 4, $intron_start + 7);
    $boundary_tree->insert($_, $intron_end - 8, $intron_end + 3);

    # cache this as it affects whether we should call something as overlapping an exon
    $tr->{_variation_effect_feature_cache}->{_has_frameshift_intron} = 1 if abs($intron_end - $intron_start) <= 12;
  }

  $tr->{_variation_effect_feature_cache}->{_intron_interval_tree} = $intron_tree;
  $tr->{_variation_effect_feature_cache}->{_intron_boundary_interval_tree} = $boundary_tree;
}

sub _exon_interval_tree {
  my $self = shift;

  my $tr = $self->transcript;

  if(!exists($tr->{_variation_effect_feature_cache}->{_exon_interval_tree})) {

    if(!$self->_can_use_interval_tree) {
      $tr->{_variation_effect_feature_cache}->{_exon_interval_tree} = undef;
    }

    else {
      my $tree = Set::IntervalTree->new();
      $tree->insert($_, $_->start - 1, $_->end) for @{$self->_exons};
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
sub _three_prime_utr {
    my $self = shift;
    
    my $tran = $self->transcript;
    
    unless(exists($tran->{_variation_effect_feature_cache}->{three_prime_utr})) {
        
        # transfer to feature slice so we don't subseq whole chromosome
        my $transferred = $tran->transfer($tran->feature_Slice());
        $tran->{_variation_effect_feature_cache}->{three_prime_utr} = $transferred->three_prime_utr();
    }
    
    return $tran->{_variation_effect_feature_cache}->{three_prime_utr};
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