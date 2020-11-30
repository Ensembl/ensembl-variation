=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2020] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY kind, either express or implied.
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

Bio::EnsEMBL::Variation::TranscriptVariationAllele

=head1 SYNOPSIS

    use Bio::EnsEMBL::Variation::TranscriptVariationAllele;
    
    my $tva = Bio::EnsEMBL::Variation::TranscriptVariationAllele->new(
        -transcript_variation   => $tv,
        -variation_feature_seq  => 'A',
        -is_reference           => 0,
    );

    print "sequence with respect to the transcript: ", $tva->feature_seq, "\n";
    print "sequence with respect to the variation feature: ", $tva->variation_feature_seq, "\n";
    print "consequence SO terms: ", (join ",", map { $_->SO_term } @{ $tva->get_all_OverlapConsequences }), "\n";
    print "amino acid change: ", $tva->pep_allele_string, "\n";
    print "resulting codon: ", $tva->codon, "\n";
    print "reference codon: ", $tva->transcript_variation->get_reference_TranscriptVariationAllele->codon, "\n";
    print "PolyPhen prediction: ", $tva->polyphen_prediction, "\n";
    print "SIFT prediction: ", $tva->sift_prediction, "\n";

=head1 DESCRIPTION

A TranscriptVariationAllele object represents a single allele of a TranscriptVariation.
It provides methods that are specific to the sequence of the allele, such as codon,
peptide etc. Methods that depend only on position (e.g. CDS start) will be found in 
the associated TranscriptVariation. Ordinarily you will not create these objects 
yourself, but instead you would create a TranscriptVariation object which will then 
construct TranscriptVariationAlleles based on the allele string of the associated
VariationFeature. 

Note that any methods that are not specific to Transcripts will be found in the 
VariationFeatureOverlapAllele superclass.

=cut

package Bio::EnsEMBL::Variation::TranscriptVariationAllele;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix qw($AA_LOOKUP);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Variation::Utils::Sequence qw(hgvs_variant_notation format_hgvs_string get_3prime_seq_offset);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(overlap within_cds within_intron stop_lost start_lost frameshift stop_retained);

use base qw(Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele Bio::EnsEMBL::Variation::BaseTranscriptVariationAllele);


our $DEBUG = 0;
our $NO_TRANSFER = 0;

sub new_fast {
    my ($class, $hashref, $strong) = @_;
    
    # swap a transcript_variation argument for a variation_feature_overlap one
    if ($hashref->{transcript_variation}) {
        $hashref->{variation_feature_overlap} = delete $hashref->{transcript_variation};
    }

    # and call the superclass

    return $class->SUPER::new_fast($hashref, $strong);
}

=head2 _return_3prime

  Description: Shifts an insertion/deletion as far as possible in the 3' direction, relative
               to the transcript.
               Called in TranscriptVariation.pm for original mapping and in 
               TranscriptVariationAllele.pm if required for HGVS
  Exceptions : none

=cut

sub _return_3prime {
    
  ## Will create a "shift_hash", containing info on precisely how the variant should be shifted when required
  my $self = shift;
  my $hgvs_only = shift;
  
  ## Return if we have already calculated a shifting hash for this allele
  return if defined($self->{shift_hash}) || ($self->is_reference && !$hgvs_only);
  
  my $tv = $self->transcript_variation;
  my $vf = $tv->base_variation_feature;
  my $var_class = $vf->var_class;
  
  ## Don't even attempt shifting if it's not an insertion or deletion
  return unless ($var_class eq 'insertion' || $var_class eq 'deletion' );

  ## If the associated database adaptor has switched HGVS shifting off, don't shift anything
  return if (defined($vf->adaptor) && defined($vf->adaptor->db)) && ($vf->adaptor->db->shift_hgvs_variants_3prime() == 0);
  
  my $vf_allele_string = $vf->allele_string;
  return if $vf_allele_string =~ /INS|DEL|HGMD|COSM|CNV/i;

  my $tr = $tv->transcript; 
  $self->initialise_unshifted_values;
  
  ## Checks to see if a simliar shift hash already exists for this variation feature.
  ## If not, a genomic shift is performed
  my $hash = $self->check_tva_shifting_hashes;
  $self->{shift_hash} = $hash if(defined($hash));
  
  ## Performs a shift in either the 5' or 3' direction depending on the strand of the transcript
  if (!defined($self->{shift_hash}))
  {
    $self->_genomic_shift($tr->strand);
    $hash = $self->check_tva_shifting_hashes;
    $self->{shift_hash} = $hash if(defined($hash))
  }

  return unless (defined($tv->cdna_start) && defined($tv->cdna_end) && defined($tv->cds_start) && defined($tv->cds_end)) 
  || (defined($tv->cdna_start_unshifted) && defined($tv->cdna_end_unshifted) && defined($tv->cds_start_unshifted) && defined($tv->cds_end_unshifted));
  
  return if ($tr->stable_id =~ /ENS/);
  my @attribs = @{$tr->get_all_Attributes()};
  
  ## Checks to see if the underlying sequence has been edited
  my @edit_attrs = grep {$_->code =~ /^_rna_edit/ && !$self->is_polyA($_)} @attribs;
  
  ## Return unless the RefSeq transcript has edited sequence - we need to check flanking sequences
  return unless scalar(@edit_attrs);
  
  my $hgvs_notation;

  ## Get all shift hashes generated for the associated variation feature
  ## If the flanking sequences are the same for (shift_length + 1) bases, we can just copy
  ## that shift hash rather than calculating another one
  my @preshifted_hashes = @{$self->base_variation_feature->{tva_shift_hashes}};
  if(scalar(@preshifted_hashes))
  {
    my ($slice_start2, $slice_end2, $slice ) = $self->_var2transcript_slice_coords($tr, $tv, $vf);

    if(defined($slice))
    {
      ## Sort through each shift hash with a corresponding transcript slice 
      ## Get the flanking sequences from the slice object and compare to previous hashes
      foreach my $shifted_obj (@preshifted_hashes)
      { 
        my $start = $var_class eq 'insertion' ? $tv->cdna_end_unshifted : $tv->cdna_start_unshifted;
        my $end = $var_class eq 'insertion' ? $tv->cdna_start_unshifted : $tv->cdna_end_unshifted;

        my $seq_to_search = defined($tr->{_variation_effect_feature_cache}->{spliced_seq}) ? $tr->{_variation_effect_feature_cache}->{spliced_seq} : $tr->seq->seq;
        my $whole_seq = substr($seq_to_search, $start - $shifted_obj->{shift_length} - 2, $end - $start + 1 + (2 * ($shifted_obj->{shift_length} + 1)));
        reverse_comp(\$whole_seq) if $tr->strand != 1;
        
        if(($shifted_obj->{type} eq 'ins' && (length($whole_seq) != ((2 * ($shifted_obj->{shift_length} + 1))))) || 
          ($shifted_obj->{type} eq 'del' && (length($whole_seq) != ((2 * ($shifted_obj->{shift_length} + 1)) + length($shifted_obj->{allele_string})))) )
        {
          ## This happens when an insertion/deletion gets too close to the transcript boundary and shifting in either direction at transcript level becomes tricky.
          $self->{shift_hash} = $shifted_obj;
          return;
        }
        
        my $pre_substr = substr($whole_seq, 0 , 1 + $shifted_obj->{shift_length});
        my $post_substr = substr($whole_seq, , -1 - $shifted_obj->{shift_length});
        if ($pre_substr eq $shifted_obj->{five_prime_flanking_seq} && $post_substr eq $shifted_obj->{three_prime_flanking_seq})
        {
          ## If both flanking sequences match, just copy over the shift hash
          $self->{shift_hash} = $shifted_obj;
          return;
        }
        
      }
    }
    else { ## This vf does not lie within that feature
      $self->{shift_hash} = $vf->{shift_hash} if (defined($vf->{shift_hash}) && $tr->strand == 1);
      $self->{shift_hash} = $vf->{shift_hash_reverse} if (defined($vf->{shift_hash_reverse}) && $tr->strand == -1);

      return;
    }
  }
  
  ## If we have no matching shift hashes already calculated, then we'll have to generate it
  
  my $edited_transcript_seq = defined($tr->{_variation_effect_feature_cache}->{spliced_seq}) ? $tr->{_variation_effect_feature_cache}->{spliced_seq} : $tr->seq->seq;
  my $start = $var_class eq 'insertion' ? $tv->cdna_end_unshifted : $tv->cdna_start_unshifted;
  my $end = $var_class eq 'insertion' ? $tv->cdna_start_unshifted : $tv->cdna_end_unshifted;

  ## To prevent having to work with too large a sequence, we get transcript sequence, get sequence +/- 1000bp from variant location, and check how far it can shift.
  ## Could be tidied up into a 'Shift 50 bases and if it's longer than that then take a larger slice'
  my $area_to_search = 1000;
  
  ## Gets flanking sequences around the variant to test if shifting is possible
  my $search_start = ($start - $area_to_search - 1) < 0 ? 0 : ($start - $area_to_search - 1);
  my $search_end = ($end + $area_to_search > length($edited_transcript_seq)) ? length($edited_transcript_seq) : ($end + $area_to_search);
  my $seqs = substr($edited_transcript_seq, $search_start, $search_end - $search_start);
  my $pre_seq = substr($seqs, 0, $start - $search_start - 1); 
  my $post_seq = substr($seqs, 0 - ($search_end - $end));
  
  my $unshifted_allele_string = $self->allele_string;
  my @allele_string = split('/', $unshifted_allele_string);
  my $hgvs_output_string = $allele_string[1];
  
  
  ## isolate correct sequence to attempt to shift
  my $seq_to_check;
  my $type;
  if ($var_class eq 'deletion')
  {
    $seq_to_check = $allele_string[0];
    $type = 'del';
  }
  elsif ($var_class eq 'insertion') 
  {
    $seq_to_check = $allele_string[1];
    $type = 'ins';
  }
  
  my $shift_length;
  my $strand = $tr->strand;
  
  ## Actually performs the shift, and provides raw data in order to create shifting hash
  ($shift_length, $seq_to_check, $hgvs_output_string, $start, $end) = $self->perform_shift($seq_to_check, $post_seq, $pre_seq, $start, $end, $hgvs_output_string, (-1 * ($strand -1))/2); 
  
  ## Creates shift_hash to attach to VF and TVA objects for 
  $self->create_shift_hash($seq_to_check, $post_seq, $pre_seq, $start, $end, $hgvs_output_string, $type, $shift_length, $strand, 0);

  return;
}

=head2 check_tva_shifting_hashes

  Description: Checks to see if a matching shifting hash already exists for another
               tva object attached to the same variation feature
  Returntype : hash
  Status     : At Risk

=cut

sub check_tva_shifting_hashes {
  my $self = shift;
  my $preshifted_hashs = $self->base_variation_feature->{tva_shift_hashes};  
  my @shift_hashes = grep {$_->{unshifted_allele_string} eq $self->allele_string  && $self->transcript->strand eq $_->{strand}} @{$preshifted_hashs} if defined($preshifted_hashs);
  if (scalar(@shift_hashes) == 1)
  {
    return $shift_hashes[0];
  }
  return undef;
}

=head2 perform_shift

  Description: Calculates the maximum possible shift in the 3' direction, provides the 
               shift length, shifted allele string, and new start/end coordinates
               Requires flanking sequences, allele string, strand and feature location
  Returntype : 5 scalars - shift length, shifted allele string, hgvs allele string, start position, end position
  Status     : At Risk

=cut

sub perform_shift {
  ## Performs the shifting calculation
  my ($self, $seq_to_check, $post_seq, $pre_seq, $var_start, $var_end, $hgvs_output_string, $reverse) = @_;
  ## get length of pattern to check 
  my $indel_length = (length $seq_to_check);
  my $shift_length = 0;
  
  ## Sets up the for loop, ensuring that the correct bases are compared depending on the strand
  my $loop_limiter = $reverse ? (length($pre_seq) - $indel_length) + 1 : (length($post_seq) - $indel_length);
  $loop_limiter = length($post_seq) if $loop_limiter < 0;
  
  for (my $n = $reverse; $n <= $loop_limiter; $n++ ){
    ## check each position in deletion/ following seq for match
    my $check_next_del;
    my $check_next_pre;
    my $hgvs_next_del;
    
    if ($reverse)
    {
      $check_next_del  = substr($seq_to_check, length($seq_to_check) -1, 1);
      $check_next_pre = substr($pre_seq, length($pre_seq) - $n, 1);
      $hgvs_next_del  = substr($hgvs_output_string, length($hgvs_output_string) -1, 1);
    }
    else{
      $check_next_del  = substr($seq_to_check, 0, 1);
      $check_next_pre = substr($post_seq, $n, 1);
      $hgvs_next_del  = substr($hgvs_output_string, 0, 1);  
    }
    if($check_next_del eq $check_next_pre){
      
      ## move position of deletion along
      $shift_length++;
      
      ## Reforms the sequences to check and the HGVS output strings for the next iteration
      if($reverse)
      {
        $seq_to_check = substr($seq_to_check, 0, length($seq_to_check) -1);
        $hgvs_output_string = substr($hgvs_output_string, 0, length($hgvs_output_string) -1);
        $seq_to_check = $check_next_del . $seq_to_check;  
        $hgvs_output_string = $hgvs_next_del . $hgvs_output_string;
      }
      else
      {
        $seq_to_check = substr($seq_to_check,1);
        $hgvs_output_string = substr($hgvs_output_string,1);
        $seq_to_check .= $check_next_del;
        $hgvs_output_string .= $hgvs_next_del;
        $var_start++;
        $var_end++;
      }
    }
    else{
      last;	    
    }
  }
  
  return $shift_length, $seq_to_check, $hgvs_output_string, $var_start, $var_end;
}

=head2 create_shift_hash

  Description: Generates a hash to attach to the TVA object to store shifting info.
               Can contain shifting info for both directions for genes with 
               transcripts on both strands
  Returntype : hash
  Status     : Stable

=cut


sub create_shift_hash
{
  my ($self, $seq_to_check, $post_seq, $pre_seq, $var_start, $var_end, $hgvs_output_string, $type, $shift_length, $strand, $genomic) = @_;
  my $vf = $self->variation_feature;
  my $five_prime_flanking_seq = substr($pre_seq, -1 - $shift_length);
  my $three_prime_flanking_seq = substr($post_seq, 0, $shift_length + 1);
  
  my @allele_string = split('/', $self->allele_string);
  reverse_comp(\$seq_to_check) if $vf->strand <0; 
  my %shift_hash = (
    "genomic" => $genomic,
    "strand" => $strand,
    "shifted_allele_string"  => $seq_to_check,
    "unshifted_allele_string" => $self->allele_string, 
    "shift_length"  => $shift_length,
    "start" => $var_start,
    "end" => $var_end,
    "type" => $type,
    "unshifted_start" => $vf->seq_region_start,
    "unshifted_end" => $vf->seq_region_end,
    "hgvs_allele_string" => $hgvs_output_string,
    "ref_orig_allele_string" => $allele_string[0],
    "alt_orig_allele_string" => $allele_string[1],
    "_hgvs_offset" => $shift_length,
    "five_prime_flanking_seq" => $five_prime_flanking_seq,
    "three_prime_flanking_seq" => $three_prime_flanking_seq,
    "allele_string" => $self->allele_string, 
  );
    $self->{shift_hash} = \%shift_hash unless $genomic;
    $vf->{shift_hash} = \%shift_hash if ($strand == 1 && $genomic);
    $vf->{shift_hash_reverse} = \%shift_hash if ($strand == -1 && $genomic);
    $vf->{tva_shift_hashes} = [] unless defined($vf->{tva_shift_hashes});
    push @{$vf->{tva_shift_hashes}}, \%shift_hash;
    
    
    return %shift_hash;
}

=head2 _genomic_shift

  Description: Performs the shifting at genomic level - all that's required for EnsEMBL 
               transcripts. Shifts 3' on either strand
  Status     : Stable

=cut

sub _genomic_shift
{
  ## Does the initial shift at the genomic level, can be on either strand
  my $self = shift;
  my $strand = shift;
  my $tv = $self->transcript_variation;
  my $vf = $tv->variation_feature;
  my $var_class = $vf->var_class;

  return unless ($var_class eq 'insertion' || $var_class eq 'deletion' );
  
  my $slice_to_shrink = $vf->slice;
  my ($slice_start, $slice_end, $var_start, $var_end) = ($slice_to_shrink->start, $slice_to_shrink->end, $vf->seq_region_start, $vf->seq_region_end );
  
  ## Gets chr slice, gets sequence +/- 1000bp from variant location, and checks how far it can shift.
  ## Could be tidied up into a 'Shift 50 bases and if it's longer than that then take a larger slice'
  my $area_to_search = 1000;
  my $orig_start = $var_start;
  my $orig_end = $var_end;
  
  my $new_slice = $slice_to_shrink->expand(0 - ($var_start - $slice_start - $area_to_search), 0 - ($slice_end - $var_end - $area_to_search));
  $new_slice = $new_slice->constrain_to_seq_region();
  
  ## Gets flanking sequences around the variant to test if shifting is possible
  my $seqs = $new_slice->seq;
  my $pre_seq = substr($seqs, 0, $area_to_search); 
  my $post_seq = substr($seqs, 0 - $area_to_search);
  
  my $unshifted_allele_string = $self->allele_string;
  my @allele_string = split('/', $unshifted_allele_string);
  $allele_string[1] = '-' if (($self->is_reference) && (scalar(@allele_string) == 1));
  my $hgvs_output_string = $allele_string[1];
  
  ## isolate correct sequence to attempt to shift
  my $seq_to_check;
  my $type;
  if ($var_class eq 'deletion')
  {
    $seq_to_check = $allele_string[0];
    $type = 'del';
  }
  elsif ($var_class eq 'insertion') 
  {
    $seq_to_check = $allele_string[1];
    $type = 'ins';
  }
  
  my $shift_length;
  reverse_comp(\$seq_to_check) if $self->variation_feature->strand <0; 
  
  ## Actually performs the shift, and provides raw data in order to create shifting hash
  ($shift_length, $seq_to_check, $hgvs_output_string, $var_start, $var_end) = $self->perform_shift($seq_to_check, $post_seq, $pre_seq, $var_start, $var_end, $hgvs_output_string, (-1 * ($strand -1))/2); 
  
  ## Creates shift_hash to attach to VF and TVA objects for 
  $self->create_shift_hash($seq_to_check, $post_seq, $pre_seq, $var_start, $var_end, $hgvs_output_string, $type, $shift_length, $strand, 1);
}

=head2 look_for_slice_start

  Description: Finds slice_start value if it can't be found in the _return_3prime method
  Status     : At Risk

=cut

sub look_for_slice_start {
  my $self = shift;
  
  my $tv = $self->transcript_variation;
  my $vf = $tv->base_variation_feature;
  my $tr = $tv->transcript;
  
  my ($slice_start, $slice_end, $slice ) = $self->_var2transcript_slice_coords($tr, $tv, $vf);
  ## set new HGVS string
  if($slice_start)
  {
    ## add cache of seq/ pos required by c, n and p 
    $self->{_slice_start} = $slice_start;
    $self->{_slice_end}   = $slice_end;
    $self->{_slice}       = $slice;
  }
  return; 
}

=head2 clear_shifting_variables

  Description: Clears shifting variables to allow for the correct values to be 
               used by the OutputFactory and future analysis (i.e. Plugins)
  Status     : At Risk

=cut

sub clear_shifting_variables {
  
  my $self = shift;
  
  my $tv = $self->transcript_variation;
  my $tr ||= $tv->transcript;

  delete($tv->{cds_end});
  delete($tv->{cds_start});
  delete($tv->{_cds_coords});
  delete($tv->{translation_start});
  delete($tv->{translation_end});
  delete($tv->{_translation_coords});
  delete($tv->{cdna_start});
  delete($tv->{cdna_end});
  delete($tv->{_cdna_coords});
  
  if(defined($self->{shift_hash}))
  {
    $tv->cds_start(undef, $tr->strand * $self->{shift_hash}->{shift_length});
    $tv->cds_end(undef, $tr->strand * $self->{shift_hash}->{shift_length});
    $tv->cdna_start(undef, $tr->strand * $self->{shift_hash}->{shift_length});
    $tv->cdna_end(undef, $tr->strand * $self->{shift_hash}->{shift_length});
    $tv->translation_start(undef, $tr->strand * $self->{shift_hash}->{shift_length});
    $tv->translation_end(undef, $tr->strand * $self->{shift_hash}->{shift_length});
  }
}


=head2 initialise_unshifted_values

  Description: Populates the unshifted values at the beginning before any shifting
               is done
  Status     : Stable

=cut

sub initialise_unshifted_values {
  
  my $self = shift;
  
  my $tv = $self->transcript_variation;
  
  unless(defined($tv->{cds_end_unshifted}))
  {  
    $tv->cds_start_unshifted();
    $tv->cds_end_unshifted();
  }
  unless(defined($tv->{cdna_end_unshifted}))
  {  
    $tv->cdna_start_unshifted();
    $tv->cdna_end_unshifted();
  }
  unless(defined($tv->{translation_end_unshifted}))
  {  
    $tv->translation_start_unshifted();
    $tv->translation_end_unshifted();
  }
}


=head2 transcript_variation

  Description: Get/set the associated TranscriptVariation
  Returntype : Bio::EnsEMBL::Variation::TranscriptVariation
  Exceptions : throws if the argument is the wrong type
  Status     : Stable

=cut

sub transcript_variation {
    my ($self, $tv) = @_;
    assert_ref($tv, 'Bio::EnsEMBL::Variation::TranscriptVariation') if $Bio::EnsEMBL::Utils::Scalar::ASSERTIONS && $tv;
    return $self->variation_feature_overlap($tv);
}

=head2 variation_feature

  Description: Get the associated VariationFeature
  Returntype : Bio::EnsEMBL::Variation::VariationFeature
  Exceptions : none
  Status     : Stable

=cut

sub variation_feature {
    my $self = shift;
    return $self->transcript_variation->variation_feature;
}

=head2 affects_peptide

  Description: Check if this changes the resultant peptide sequence
  Returntype : boolean
  Exceptions : None
  Caller     : general
  Status     : At Risk

=cut

sub affects_peptide {
  my $self = shift;
  return scalar grep { $_->SO_term =~ /stop|missense|frameshift|inframe|initiator/ } @{$self->get_all_OverlapConsequences}; 
}

=head2 pep_allele_string

  Description: Return a '/' delimited string of the reference peptide and the 
               peptide resulting from this allele, or a single peptide if this
               allele does not change the peptide (e.g. because it is synonymous)
  Returntype : string or undef if this allele is not in the CDS
  Exceptions : none
  Status     : Stable

=cut

sub pep_allele_string {
    my ($self) = @_;

    my $pep = $self->peptide;
    
    return undef unless $pep;
    
    my $ref_pep = $self->transcript_variation->get_reference_TranscriptVariationAllele->peptide;

    return undef unless $ref_pep;
    
    return $ref_pep ne $pep ? $ref_pep.'/'.$pep : $pep;
}

=head2 codon_allele_string

  Description: Return a '/' delimited string of the reference codon and the 
               codon resulting from this allele 
  Returntype : string or undef if this allele is not in the CDS
  Exceptions : none
  Status     : Stable

=cut

sub codon_allele_string {
    my ($self) = @_;
    
    my $codon = $self->codon;
    
    return undef unless $codon;
    
    my $ref_codon = $self->transcript_variation->get_reference_TranscriptVariationAllele->codon;
    
    return $ref_codon.'/'.$codon;
}

=head2 display_codon_allele_string

  Description: Return a '/' delimited string of the reference display_codon and the 
               display_codon resulting from this allele. The display_codon identifies
               the nucleotides affected by this variant in UPPER CASE and other 
               nucleotides in lower case
  Returntype : string or undef if this allele is not in the CDS
  Exceptions : none
  Status     : At Risk

=cut

sub display_codon_allele_string {
    my ($self) = @_;
    
    my $display_codon = $self->display_codon;
    
    return undef unless $display_codon;
        
    my $ref_tva = $self->transcript_variation->get_reference_TranscriptVariationAllele;
    $ref_tva->{shift_hash} = $self->{shift_hash};
    
    my $ref_display_codon = $ref_tva->display_codon;
    
    return undef unless $ref_display_codon;
    
    return $ref_display_codon.'/'.$display_codon;
}

=head2 peptide

  Description: Return the amino acid sequence that this allele is predicted to result in
  Returntype : string or undef if this allele is not in the CDS or is a frameshift
  Exceptions : none
  Status     : Stable

=cut

sub peptide {
  my ($self, $peptide) = @_;

  $self->{peptide} = $peptide if $peptide;

  unless(exists($self->{peptide})) {

    $self->{peptide} = undef;

    return $self->{peptide} unless $self->seq_is_unambiguous_dna;

    if(my $codon = $self->codon) {

      # the codon method can set the peptide in some circumstances 
      # so check here before we try an (expensive) translation
      return $self->{peptide} if $self->{peptide};

      my $tv = $self->base_variation_feature_overlap;

      # for mithocondrial dna we need to to use a different codon table
      my $codon_table = $tv->_codon_table;

      # check the cache
      my $pep_cache = $main::_VEP_CACHE->{pep}->{$codon_table} ||= {};
      if(!($self->{is_reference} && scalar @{$tv->_seq_edits}) && ($self->{peptide} = $pep_cache->{$codon})) {
        return $self->{peptide};
      }
      
      # translate the codon sequence to establish the peptide allele
      
      # allow for partial codons - split the sequence into whole and partial
      # e.g. AAAGG split into AAA and GG            
      my $whole_codon   = substr($codon, 0, int(length($codon) / 3) * 3);
      my $partial_codon = substr($codon, int(length($codon) / 3) * 3);

      my $pep = '';
      
      if($whole_codon) {          
          my $codon_seq = Bio::Seq->new(
            -seq        => $whole_codon,
            -moltype    => 'dna',
            -alphabet   => 'dna',
            );
          
          $pep .= $codon_seq->translate(undef, undef, undef, $codon_table)->seq;
        }

      # apply any seq edits?
      my $have_edits = 0;

      if($self->{is_reference}) {
        my $seq_edits = $tv->_seq_edits;
        
        if(scalar @$seq_edits) {

          # get TV coords, switch if necessary
          my ($tv_start, $tv_end) = ($tv->translation_start, $tv->translation_end);
          ($tv_start, $tv_end) = ($tv_end, $tv_start) if $tv_start > $tv_end;
          
          # get all overlapping seqEdits
          SE: foreach my $se(grep {overlap($tv_start, $tv_end, $_->start, $_->end)} @$seq_edits) {
            my ($se_start, $se_end, $alt) = ($se->start, $se->end, $se->alt_seq);
            my $se_alt_seq_length = length($alt);
            $have_edits = 1;
            
            # loop over each overlapping pos
            foreach my $tv_pos(grep {overlap($_, $_, $se_start, $se_end)} ($tv_start..$tv_end)) {

              # in some cases, the sequence edit can shorten the protein
              # this means our TV can fall outside the range of the edited protein
              # therefore for safety jump out
              if($tv_pos - $se_start >= $se_alt_seq_length) {
                return $self->{peptide} = undef;
              }

              # apply edit, adjusting for string position
              substr($pep, $tv_pos - $tv_start, 1) = substr($alt, $tv_pos - $se_start, 1);
            }
          }
        }
      }
      
      if($partial_codon && $pep ne '*') {
        $pep .= 'X';
      }

      $pep ||= '-';

      $pep_cache->{$codon} = $pep if length($codon) <= 3 && !$have_edits;

      $self->{peptide} = $pep;
    }
  }

  return $self->{peptide};
}

=head2 codon

  Description: Return the codon sequence that this allele is predicted to result in
  Returntype : string or undef if this allele is not in the CDS or is a frameshift
  Exceptions : none
  Status     : Stable

=cut

sub codon {
  my ($self, $codon) = @_;
  $self->{codon} = $codon if defined $codon;

  unless(exists($self->{codon})) {

    $self->{codon} = undef;
  
    my $tv = $self->base_variation_feature_overlap;
    
    my $shifting_offset = 0;
    my $tr = $tv->transcript;
  
    $shifting_offset = (defined($self->{shift_hash}) && defined($self->{shift_hash}->{shift_length})) ? $self->{shift_hash}->{shift_length} : 0;

    my ($tv_tr_start, $tv_tr_end) = ($tv->translation_start, $tv->translation_end);

    unless($tv_tr_start && $tv_tr_end && $self->seq_is_dna) {
      return $self->{codon};
    }

    # try to calculate the codon sequence
    $self->shift_feature_seqs unless $shifting_offset == 0;
    my $seq = $self->feature_seq;
    $seq = '' if $seq eq '-';
    
    # calculate necessary coords and lengths
    
    my $codon_cds_start = $tv_tr_start * 3 - 2;
    my $codon_cds_end   = $tv_tr_end * 3;
    my $codon_len       = $codon_cds_end - $codon_cds_start + 1;
    unless($shifting_offset == 0)
    {
      delete($tv->{cds_end});
      delete($tv->{cds_start});
      delete($tv->{_cds_coords});
    }
    return $self->{codon} if !defined($tv->cds_end(undef, $tr->strand * $shifting_offset)) || !defined($tv->cds_start(undef, $tr->strand * $shifting_offset));
    my $vf_nt_len       = $tv->cds_end(undef, $tr->strand * $shifting_offset) - $tv->cds_start(undef, $tr->strand * $shifting_offset) + 1;
    my $allele_len      = $self->seq_length;
    unless($shifting_offset == 0)
    {
      delete($tv->{cds_end});
      delete($tv->{cds_start});
      delete($tv->{_cds_coords});
    }
    my $cds;
    if ($allele_len != $vf_nt_len) {
      if (abs($allele_len - $vf_nt_len) % 3) {
        # this is a frameshift variation, we don't attempt to 
        # calculate the resulting codon or peptide change as this 
        # could get quite complicated 
        # return undef;
      }

      ## Bioperl Seq object
      my $cds_obj = $self->_get_alternate_cds();
      return undef unless defined($cds_obj);
      $cds = $cds_obj->seq();
    }

    else {
      # splice the allele sequence into the CDS
      $cds = $tv->_translateable_seq;

      substr($cds, $tv->cds_start(undef, $tr->strand * $shifting_offset) -1, $vf_nt_len) = $seq;
    }

    # and extract the codon sequence
    my $codon = substr($cds, $codon_cds_start-1, $codon_len + ($allele_len - $vf_nt_len));
    
    if (length($codon) < 1) {
      $self->{codon}   = '-';
      $self->{peptide} = '-';
    }
    else {
       $self->{codon} = $codon;
    }
  }
  
  return $self->{codon};
}

=head2 display_codon

  Description: Return the codon sequence that this allele is predicted to result in
               with the affected nucleotides identified in UPPER CASE and other 
               nucleotides in lower case
  Returntype : string or undef if this allele is not in the CDS or is a frameshift
  Exceptions : none
  Status     : Stable

=cut

sub display_codon {
  my $self = shift;

  unless(exists($self->{_display_codon})) {

    # initialise so it doesn't get called again
    $self->{_display_codon} = undef;

    if(my $codon = $self->codon) {

      my $display_codon = lc $self->codon;

      if(my $codon_pos = $self->transcript_variation->codon_position) {

        # if this allele is an indel then just return all lowercase
        if ($self->feature_seq ne '-') {
            
          # codon_position is 1-based, while substr assumes the string starts at 0          
          my $pos = $codon_pos - 1;

          my $len = length $self->feature_seq;

          substr($display_codon, $pos, $len) = uc substr($display_codon, $pos, $len);
        }
      }

      $self->{_display_codon} = $display_codon;
    }
  }

  return $self->{_display_codon};
}

=head2 polyphen_prediction

  Description: Return the qualitative PolyPhen-2 prediction for the effect of this allele.
               (Note that we currently only have PolyPhen predictions for variants that 
               result in single amino acid substitutions in human)
  Returntype : string (one of 'probably damaging', 'possibly damaging', 'benign', 'unknown')
               if this is a missense change and a prediction is available, undef
               otherwise
  Exceptions : none
  Status     : At Risk

=cut

sub polyphen_prediction {
    my ($self, $classifier, $polyphen_prediction) = @_;
    
    $classifier ||= 'humvar';
    
    my $analysis = "polyphen_${classifier}";
    
    $self->{$analysis}->{prediction} = $polyphen_prediction if $polyphen_prediction;
    
    unless (defined $self->{$analysis}->{prediction}) {
        my ($prediction, $score) = $self->_protein_function_prediction($analysis);
        $self->{$analysis}->{score} = $score;
        $self->{$analysis}->{prediction} = $prediction;
    }
    
    return $self->{$analysis}->{prediction};
}

=head2 polyphen_score

  Description: Return the PolyPhen-2 probability that this allele is deleterious (Note that we 
               currently only have PolyPhen predictions for variants that result in single 
               amino acid substitutions in human)
  Returntype : float between 0 and 1 if this is a missense change and a prediction is 
               available, undef otherwise
  Exceptions : none
  Status     : At Risk

=cut

sub polyphen_score {
    my ($self, $classifier, $polyphen_score) = @_;
    
    $classifier ||= 'humvar';

    my $analysis = "polyphen_${classifier}";
    
    $self->{$analysis}->{score} = $polyphen_score if defined $polyphen_score;

    unless (defined $self->{$analysis}->{score}) {
        my ($prediction, $score) = $self->_protein_function_prediction($analysis);
        $self->{$analysis}->{score} = $score;
        $self->{$analysis}->{prediction} = $prediction;
    }

    return $self->{$analysis}->{score};
}

=head2 sift_prediction

  Description: Return the qualitative SIFT prediction for the effect of this allele.
               (Note that we currently only have SIFT predictions for variants that 
               result in single amino acid substitutions in human)
  Returntype : string (one of 'tolerated', 'deleterious') if this is a missense 
               change and a prediction is available, undef otherwise
  Exceptions : none
  Status     : At Risk

=cut

sub sift_prediction {
    my ($self, $sift_prediction) = @_;
    
    $self->{sift_prediction} = $sift_prediction if $sift_prediction;
    
    unless (defined $self->{sift_prediction}) {
        my ($prediction, $score) = $self->_protein_function_prediction('sift');
        $self->{sift_score} = $score;
        $self->{sift_prediction} = $prediction unless $self->{sift_prediction};
    }
    
    return $self->{sift_prediction};
}

=head2 sift_score

  Description: Return the SIFT score for this allele (Note that we currently only have SIFT 
               predictions for variants that result in single amino acid substitutions in human)
  Returntype : float between 0 and 1 if this is a missense change and a prediction is 
               available, undef otherwise
  Exceptions : none
  Status     : At Risk

=cut

sub sift_score {
    my ($self, $sift_score) = @_;

    $self->{sift_score} = $sift_score if defined $sift_score;

    unless (defined $self->{sift_score}) {
        my ($prediction, $score) = $self->_protein_function_prediction('sift');
        $self->{sift_score} = $score;
        $self->{sift_prediction} = $prediction;
    }

    return $self->{sift_score};
}

=head2 cadd_prediction

  Description: Return the qualitative CADD prediction for the effect of this allele.
               (Note that we currently only have predictions for variants that 
               result in single amino acid substitutions in human)
  Returntype : string (one of 'likely benign', 'likely deleterious') if this is a missense 
               change and a prediction is available, undef otherwise. Predictions
               are assigned based on CADD PHRED scores. CADD PHRED scores greater or
               equal to 15 are considered likely deleterious.  
  Exceptions : none
  Status     : At Risk

=cut

sub cadd_prediction {
  my ($self, $cadd_prediction) = @_;
  return $self->_prediction('cadd_prediction', $cadd_prediction);
}

=head2 dbnsfp_revel_prediction

  Description: Return the qualitative REVEL prediction for the effect of this allele.
               (Note that we currently only have predictions for variants that 
               result in single amino acid substitutions in human)
  Returntype : string (one of 'likely_disease_causing', 'likely_not_disease_causing')
               if this is a missense change and a prediction is available, undef otherwise.
               We chose 0.5 as the threshold to assign the predictions. From the REVEL paper:
               For example, 75.4% of disease mutations but only 10.9% of neutral variants
               have a REVEL score above 0.5, corresponding to a sensitivity of 0.754 and
               specificity of 0.891. 
  Exceptions : none
  Status     : At Risk

=cut

sub dbnsfp_revel_prediction {
  my ($self, $dbnsfp_revel_prediction) = @_;
  return $self->_prediction('dbnsfp_revel_prediction', $dbnsfp_revel_prediction);
}

=head2 dbnsfp_meta_lr_prediction

  Description: Return the qualitative MetaLR prediction for the effect of this allele.
               (Note that we currently only have predictions for variants that 
               result in single amino acid substitutions in human)
  Returntype : string (one of 'tolerated', 'damaging').
               The score cutoff between "D" and "T" is 0.5.
  Exceptions : none
  Status     : At Risk

=cut

sub dbnsfp_meta_lr_prediction {
  my ($self, $dbnsfp_meta_lr_prediction) = @_;
  return $self->_prediction('dbnsfp_meta_lr_prediction', $dbnsfp_meta_lr_prediction);
}

=head2 dbnsfp_mutation_assessor_prediction

  Description: Return the qualitative MutationAssessor prediction for the effect of this allele.
               (Note that we currently only have predictions for variants that 
               result in single amino acid substitutions in human)
  Returntype : string (one of 'high', 'medium', 'low', 'neutral').
  Exceptions : none
  Status     : At Risk

=cut

sub dbnsfp_mutation_assessor_prediction {
  my ($self, $dbnsfp_mutation_assessor_prediction) = @_;
  return $self->_prediction('dbnsfp_mutation_assessor_prediction', $dbnsfp_mutation_assessor_prediction);
}

=head2 _prediction

  Description: Return prediction for specified score type.
  Returntype : float 
  Exceptions : none
  Status     : At Risk

=cut

sub _prediction {
  my ($self, $prediction_type, $prediction) = @_;
  $self->{$prediction_type} = $prediction if $prediction;

  unless (defined $self->{$prediction_type}) {
    my $analysis = $prediction_type;
    $analysis =~ s/_prediction//;
    my ($prediction, $score) = $self->_protein_function_prediction($analysis);
    my $score_type = $prediction_type;
    $score_type =~ s/_prediction/_score/;
    $self->{$score_type} = $score;
    $self->{$prediction_type} = $prediction;
  }
  return $self->{$prediction_type};
}

=head2 cadd_score

  Description: Return the CADD PHRED score for this allele.
  Returntype : float if this is a missense change and a prediction is available, undef otherwise
  Exceptions : none
  Status     : At Risk

=cut

sub cadd_score {
  my ($self, $cadd_score) = @_;
  return $self->_score('cadd_score');
}

=head2 dbnsfp_revel_score

  Description: Return the REVEL score for this allele. The score is retrieved from dbNSFP. (We only
               have predictions for variants that result in single amino acid substitutions in human)
  Returntype : float between 0 and 1 if this is a missense change and a prediction is 
               available, undef otherwise
  Exceptions : none
  Status     : At Risk

=cut

sub dbnsfp_revel_score {
  my ($self, $dbnsfp_revel_score) = @_;
  return $self->_score('dbnsfp_revel_score');
}

=head2 dbnsfp_meta_lr_score

  Description: Return the MetaLR score for this allele. The score is retrieved from dbNSFP. (We only
               have predictions for variants that result in single amino acid substitutions in human)
  Returntype : float if this is a missense change and a prediction is 
               available, undef otherwise
  Exceptions : none
  Status     : At Risk

=cut

sub dbnsfp_meta_lr_score {
  my ($self, $dbnsfp_meta_lr_score) = @_;
  return $self->_score('dbnsfp_meta_lr_score', $dbnsfp_meta_lr_score);
}

=head2 dbnsfp_mutation_assessor_score

  Description: Return the MutationAssessor score for this allele. The score is retrieved from dbNSFP. (We only
               have predictions for variants that result in single amino acid substitutions in human)
  Returntype : float if this is a missense change and a prediction is 
               available, undef otherwise
  Exceptions : none
  Status     : At Risk

=cut

sub dbnsfp_mutation_assessor_score {
  my ($self, $dbnsfp_mutation_assessor_score) = @_;
  return $self->_score('dbnsfp_mutation_assessor_score', $dbnsfp_mutation_assessor_score);
}


=head2 _score

  Description: Return score for specified score type.
  Returntype : float 
  Exceptions : none
  Status     : At Risk

=cut

sub _score {
  my ($self, $score_type, $score) = @_; 
  $self->{$score_type} = $score if defined $score;

  unless (defined $self->{$score_type}) {
      my $analysis = $score_type;
      $analysis =~ s/_score//;
      my ($prediction, $score) = $self->_protein_function_prediction($analysis);
      my $prediction_type = $score_type;
      $prediction_type =~ s/_score/_prediction/;
      $self->{$score_type} = $score;
      $self->{$prediction_type} = $prediction;
  }
  return $self->{$score_type};
}

sub _protein_function_prediction {
    my ($self, $analysis) = @_;

    # we can only get results for variants that cause a single amino acid substitution, 
    # so check the peptide allele string first

    if ($self->pep_allele_string && $self->pep_allele_string =~ /^[A-Z]\/[A-Z]$/ && defined $AA_LOOKUP->{$self->peptide}) {
        if (my $matrix = $self->transcript_variation->_protein_function_predictions($analysis)) {
          
            # temporary fix - check $matrix is not an empty hashref
            if(ref($matrix) && ref($matrix) eq 'Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix') {
            
                my ($prediction, $score) = $matrix->get_prediction(
                    $self->transcript_variation->translation_start,
                    $self->peptide,
                );

                return wantarray ? ($prediction, $score) : $prediction;
            }
        }
    }
    
    return undef;
}

=head2 shift_feature_seqs

  Description: Updates the feature_seq and variation_feature_seq values associated with this TVA to be shifted
  Returntype : string 
  Exceptions : none
  Status     : At Risk

=cut

sub shift_feature_seqs {
    my $self = shift;
    if(defined($self->{shift_hash}) && !$self->{shifted_feature_seqs})
    {
      my $vf_seq = $self->{variation_feature_seq};
      my $f_seq = $self->{feature_seq};
      
      my $shift_length = $self->{shift_hash}->{shift_length} ||= 0;
      $shift_length = $self->seq_length - $shift_length if $self->transcript->strand == -1;
      for (my $n = 0; $n < $shift_length; $n++ ){
        ## check each position in deletion/ following seq for match
        $vf_seq = substr($vf_seq, 1) . substr($vf_seq, 0, 1) if defined($vf_seq);
        $f_seq = substr($f_seq, 1) . substr($f_seq, 0, 1) if defined($f_seq);
      } 
      $self->variation_feature_seq($vf_seq);
      $self->feature_seq($f_seq);
      $self->{shifted_feature_seqs} = 1;
    }   
}


=head2 hgvs_genomic

  Description: Return a string representing the genomic-level effect of this allele in HGVS format
  Returntype : string 
  Exceptions : none
  Status     : At Risk

=cut

sub hgvs_genomic {
    return _hgvs_generic(@_,'genomic');
}


=head2 hgvs_transcript

  Description: Return a string representing the CDS-level effect of this allele in HGVS format
  Returntype : string or undef if this allele is not in the CDS
  Exceptions : none
  Status     : At Risk

=cut
    
    
sub hgvs_transcript {
  my $self = shift;
  my $notation = shift;
  my $no_shift = shift;
  
  ##### set if string supplied
  $self->{hgvs_transcript} = $notation   if defined $notation;

  ##### return if held 
  return $self->{hgvs_transcript}        if defined $self->{hgvs_transcript};

  ## This is cached to allow long form HGVS to be created
  ## Set as undefined here to avoid re-calling if hgvs_transcript annotation not possible
  $self->{hgvs_t_ref} = undef;

  my $tv = $self->base_variation_feature_overlap;

  ### don't attempt to recalculate if field is NULL from DB
  return $self->{hgvs_transcript}        if exists $self->{hgvs_transcript} && defined($tv->dbID);


  ### don't try to handle odd characters
  return undef if $self->variation_feature_seq() =~ m/[^ACGT\-]/ig;

  ### no result for reference allele
  return undef if $self->is_reference == 1;

  ### else evaluate
  my $tr = $tv->transcript;
  my $tr_stable_id = $tr->stable_id;
  my $vf = $tv->base_variation_feature;
    
  ### get reference sequence strand
  my $refseq_strand = $tr->strand();

  my $var_name = $vf->variation_name();

  if($DEBUG ==1){    
  	print "\nHGVS transcript: Checking ";
  	print " var_name $var_name " if defined $var_name ;
  	print " refseq strand => $refseq_strand"  if defined $refseq_strand;
  	print " seq name : " . $tr_stable_id  
  	    if defined $tr_stable_id ;
  	print " var strand " . $vf->strand()  
  	    if defined $vf->strand();
  	print " vf st " . $vf->strand()   
  		if defined $vf->strand()  ;
  	print " seqname: " . $vf->seqname()  ." seq: " . $self->variation_feature_seq  
  	    if defined  $vf->seqname();
  	print "\n";
  }

  my $hgvs_notation ; ### store components of HGVS string in hash

  my $variation_feature_sequence;
  my $adaptor_shifting_flag = 1;
  ## Check previous shift_hgvs_variants_3prime flag and act accordingly
  if (defined($vf->adaptor) && defined($vf->adaptor->db)) {
    $adaptor_shifting_flag = $vf->adaptor->db->shift_hgvs_variants_3prime();
  }
  elsif(defined($Bio::EnsEMBL::Variation::DBSQL::DBAdaptor::DEFAULT_SHIFT_HGVS_VARIANTS_3PRIME)){
    $adaptor_shifting_flag = $Bio::EnsEMBL::Variation::DBSQL::DBAdaptor::DEFAULT_SHIFT_HGVS_VARIANTS_3PRIME;
  }
  
  my $hash_already_defined = defined($self->{shift_hash}); 
  ## Perform HGVS shift even if no_shift is on
  $self->_return_3prime(1) unless ($adaptor_shifting_flag == 0);

  $variation_feature_sequence = $self->variation_feature_seq();
  $variation_feature_sequence = $self->{shift_hash}->{hgvs_allele_string} if defined($self->{shift_hash}) && $vf->var_class() eq 'insertion' && ($adaptor_shifting_flag != 0);
 
  my $offset_to_add = defined($self->{shift_hash}) ? $self->{shift_hash}->{_hgvs_offset} : 0;# + ($no_shift ? 0 : (0 - $self->{_hgvs_offset}) );
  $self->{_hgvs_offset} = $offset_to_add;
  ## delete the shifting hash if we generated it for HGVS calculations
  delete($self->{shift_hash}) unless $hash_already_defined;
 
  ## return if a new transcript_variation_allele is not available - variation outside transcript
  return undef unless defined $self->base_variation_feature_overlap;
  $self->look_for_slice_start unless (defined  $self->{_slice_start});
  
  unless (defined  $self->{_slice_start} ){
  	print "Exiting hgvs_transcript: no slice start position for $var_name in trans" . $tr_stable_id . "\n " if $DEBUG == 1 ;
  	return undef;
  }
  ## this may be different to the input one for insertions/deletions
    print "vfs: $variation_feature_sequence &  $self->{_slice_start} -> $self->{_slice_end}\n" if $DEBUG ==1;
  if($variation_feature_sequence && $vf->strand() != $refseq_strand) {    
    reverse_comp(\$variation_feature_sequence) ;
  };
  ## delete consequences if we have an offset. This is only in here for when we want HGVS to shift but not consequences.
  ## TODO add no_shift flag test
  delete($self->{_predicate_cache}) if $self->transcript_variation->{shifted} && $offset_to_add != 0; 
  print "sending alt: $variation_feature_sequence &  $self->{_slice_start} -> $self->{_slice_end} for formatting\n" if $DEBUG ==1;
  
  return undef if (($self->{_slice}->end - $self->{_slice}->start + 1) < ($self->{_slice_end} + $offset_to_add));
  #return undef if (length($self->{_slice}->seq()) < ($self->{_slice_end} + $offset_to_add));
  $hgvs_notation = hgvs_variant_notation(
    $variation_feature_sequence,    ### alt_allele,
    $self->{_slice}->seq(),                             ### using this to extract ref allele
    $self->{_slice_start} + $offset_to_add,
    $self->{_slice_end} + $offset_to_add,
    "",
    "",
    $var_name 
  );
  
  ### This should not happen
  unless($hgvs_notation->{'type'}){
    #warn "Error - not continuing; no HGVS annotation\n";
    return undef;
  } 

  ## check for the same bases in ref and alt strings before or after the variant
  $hgvs_notation = _clip_alleles($hgvs_notation) unless $hgvs_notation->{'type'} eq 'dup';

  print "hgvs transcript type : " . $hgvs_notation->{'type'} . "\n" if $DEBUG == 1;    
  print "Got type: " . $hgvs_notation->{'type'} ." $hgvs_notation->{'ref'} -> $hgvs_notation->{'alt'}\n" if $DEBUG == 1;

  ### create reference name - transcript name & seq version
  my $stable_id = $tr_stable_id;    
  $stable_id .= "." . $tr->version() 
     unless ($stable_id =~ /\.\d+$/ || $stable_id =~ /LRG/); ## no version required for LRG's
  $hgvs_notation->{'ref_name'} = $stable_id;


  ### get position relative to transcript features [use HGVS coords not variation feature coords due to dups]
  # avoid doing this twice if start and end are the same
  my $same_pos = $hgvs_notation->{start} == $hgvs_notation->{end};
  
  ## Mismatches between refseq transcripts and the reference genome are tracked in transcript attributes
  my @attribs = @{$tr->get_all_Attributes()};
  my @edit_attrs = grep {$_->code =~ /^_rna_edit/} @attribs;

  my $misalignment_offset = 0;
  $misalignment_offset = $self->get_misalignment_offset(\@edit_attrs) if (scalar(@edit_attrs) && (substr($tr->stable_id, 0,3) eq 'NM_' || substr($tr->stable_id, 0,3) eq 'XM_'));
  
  if ($vf->var_class eq 'SNP' && defined($self->{pre_consequence_predicates}) && $self->{pre_consequence_predicates}->{exon} && defined($tv->cds_start) && defined($tv->cds_end)) {
    $hgvs_notation->{start} = $tv->cds_start;
    $hgvs_notation->{end}   = $hgvs_notation->{start};
  }
  else{
    $hgvs_notation->{start} = $self->_get_cDNA_position( $hgvs_notation->{start} + $misalignment_offset);
    $hgvs_notation->{end}   = $same_pos ? $hgvs_notation->{start} : $self->_get_cDNA_position( $hgvs_notation->{end} + $misalignment_offset );
  }
  return undef unless defined  $hgvs_notation->{start}  && defined  $hgvs_notation->{end} ;

  # Make sure that start is always less than end
  my ($exon_start_coord, $intron_start_offset) = $hgvs_notation->{start} =~ m/(\-?[0-9]+)\+?(\-?[0-9]+)?/;
  my ($exon_end_coord,   $intron_end_offset)   = $same_pos ? ($exon_start_coord, $intron_start_offset) : $hgvs_notation->{end} =~ m/(\-?[0-9]+)\+?(\-?[0-9]+)?/;
  $intron_start_offset ||= 0;
  $intron_end_offset   ||= 0;
  print "pre pos sort : $hgvs_notation->{start},$hgvs_notation->{end}\n" if $DEBUG ==1;
  ($hgvs_notation->{start},$hgvs_notation->{end}) = ($hgvs_notation->{end},$hgvs_notation->{start}) if (
    (($exon_start_coord  > $exon_end_coord) || ($exon_start_coord  == $exon_end_coord && $intron_start_offset > $intron_end_offset) ) &&
    $hgvs_notation->{end} !~/\*/
  );
  print "post pos sort: $hgvs_notation->{start},$hgvs_notation->{end}\n" if $DEBUG ==1;

  # save these to be able to report exon coordinates and intron distances
  $self->{_hgvs_exon_start_coordinate} = $exon_start_coord;
  $self->{_hgvs_intron_start_offset} = $intron_start_offset;
  $self->{_hgvs_exon_end_coordinate} = $exon_end_coord;
  $self->{_hgvs_intron_end_offset} = $intron_end_offset;


  if($tr->cdna_coding_start()){
    $hgvs_notation->{'numbering'} = "c";  ### set 'c' if transcript is coding 
  }    
  else{
    $hgvs_notation->{'numbering'} = "n";  ### set 'n' if transcript non-coding 
  }

  ### generic formatting 
  print "pre-format $hgvs_notation->{alt}\n" if $DEBUG ==1;
  $self->{hgvs_transcript} = format_hgvs_string( $hgvs_notation);
  if($DEBUG ==1){ print "HGVS notation: " . $self->{hgvs_transcript} . " \n"; }

  ## save the HGVS style reference sequence in case the other long form of HGVS is required
  $self->{hgvs_t_ref} = ($hgvs_notation->{type} eq 'dup') ? $hgvs_notation->{alt} :  $hgvs_notation->{ref};

  return $self->{hgvs_transcript}; 
}

=head2 is_polyA

  Description: Tests RefSeq transcript attribute to see if it's a poly-A tail
  Returntype : boolean value
  Status     : Stable

=cut

sub is_polyA {
  my $self = shift;
  my $attr = shift;
  my @spl_value = split(/ /, $attr->{value});
  if((scalar(@spl_value) eq 3) && (length($self->transcript->seq->seq) <= 1 + (length($spl_value[2]) + $spl_value[0])))
  {
      return 1;
  }
  return 0;
}

=head2 get_misalignment_offset

  Description: Calculate total offset created by RefSeq bam alignment 
  Returntype : scalar value
  Status     : Stable

=cut

sub get_misalignment_offset {
  my $self = shift;
  my $attrs = shift;
  my $mlength = 0;
  
  ## For each refseq edit transcript attribute given, isolate the insertions and deletions found 
  ## before the given variant and sum their lengths
  foreach my $attr (@{$attrs})
  {
    next if defined($attr->{description}) && $attr->description =~ /op=X/;
    ## Find out whether insertion or deletion, and get length  
    my @split_val = split(/ /,$attr->{value});
    next if (scalar(@split_val) eq 3) && ($split_val[1] - $split_val[0] + 1) eq length($split_val[2]);
    my $type = (scalar(@split_val) eq 3) ? 'ins' : 'del';
    my $var_location_start = $self->transcript_variation->cdna_start;
    $var_location_start = 0 unless defined($var_location_start);
    next if $var_location_start < $split_val[0];

    if ($type eq 'del')
    {
      $mlength += -1 - ($split_val[1] - $split_val[0]);
    }
    else{
      $mlength += length($split_val[2]);
    }
  }
   
  $self->{refseq_misalignment_offset} = $mlength;
  return $mlength;
  
}


=head2 hgvs_transcript_reference

  Description: Return a string representing the reference sequence as could be used in HGVS notation
               Useful for deletions where the recommended format formerly used the reference sequence
               but no longer does.
  Returntype : string or undef if no reference allele is available.
  Exceptions : none
  Status     : At Risk

=cut

sub hgvs_transcript_reference {

  my $self = shift;

  $self->hgvs_transcript() unless exists $self->{hgvs_t_ref};
  return $self->{hgvs_t_ref};

}


=head2 hgvs_protein

  Description: Return a string representing the protein-level effect of this allele in HGVS format
  Returntype : string or undef if this allele is not in the CDS 
  Exceptions : none
  Status     : At Risk
  
=cut

sub hgvs_protein {
  my $self     = shift;
  my $notation = shift;  
  my $hgvs_notation;

  if($DEBUG == 1){
    print "\nStarting hgvs_protein with "; 
    print " var: " . $self->transcript_variation->variation_feature->variation_name() 
      if defined $self->transcript_variation->variation_feature->variation_name() ;
    print   " trans: " . $self->transcript_variation->transcript->display_id() 
      if defined $self->transcript_variation->transcript->display_id() ;
    print "\n";
  }

  ### set if string supplied
  $self->{hgvs_protein} = $notation  if defined $notation;


  ### return if set
  return $self->{hgvs_protein}       if defined $self->{hgvs_protein} ;
  
  ### don't attempt to recalculate if field is NULL from DB
  return $self->{hgvs_protein}       if exists $self->{hgvs_protein} && defined($self->transcript_variation->dbID);
  
  ### don't try to handle odd characters
  return undef if $self->variation_feature_seq() =~ m/[^ACGT\-]/ig;

  ### no HGVS annotation for reference allele
  return undef if $self->is_reference();
  print "checking pos with hgvs prot\n" if $DEBUG ==1;

  ## HGVS requires the variant be located in as 3' position as possible

  ## return if a new transcript_variation_allele is not available - variation outside transcript
  my $tv = $self->base_variation_feature_overlap;
  return undef unless defined $tv;

  my $vf = $tv->base_variation_feature;
  my $tr          = $tv->transcript;

  my $adaptor_shifting_flag = 1;
  ## Check previous shift_hgvs_variants_3prime flag and act accordingly
  if (defined($vf->adaptor) && defined($vf->adaptor->db)) {
    $adaptor_shifting_flag = $vf->adaptor->db->shift_hgvs_variants_3prime();
  }
  elsif(defined($Bio::EnsEMBL::Variation::DBSQL::DBAdaptor::DEFAULT_SHIFT_HGVS_VARIANTS_3PRIME)){
    $adaptor_shifting_flag = $Bio::EnsEMBL::Variation::DBSQL::DBAdaptor::DEFAULT_SHIFT_HGVS_VARIANTS_3PRIME;
  }
  
  ## Check to see if the shift_hash is already defined, allowing us to remove it from associated $tva objects when we only want to shift HGVS
  my $hash_already_defined = defined($self->{shift_hash});
  ## Perform HGVS shift even if no_shift is on - only prevent shifting if shift_hgvs_variants_3prime() has been switched off.
  $self->_return_3prime(1) unless ($adaptor_shifting_flag == 0);


  my $pre         = $self->_pre_consequence_predicates;
  my $shifting_offset = 0;
  if($tr->strand() > 0) {
    $shifting_offset = (defined($self->{shift_hash}) && defined($self->{shift_hash}->{shift_length})) ? $self->{shift_hash}->{shift_length} : 0;
  }
  elsif($tr->strand < 0) {
    $shifting_offset = (defined($self->{shift_hash}) && defined($self->{shift_hash}->{shift_length})) ? 0 - $self->{shift_hash}->{shift_length} : 0;
  }

  ### no HGVS protein annotation for variants outside translated region 
  if(defined($self->{shift_hash}) && defined($self->{shift_hash}->{shift_length}) && $self->{shift_hash}->{shift_length} != 0)
  {
    delete($tv->{translation_start});
    delete($tv->{translation_end});
    delete($tv->{_translation_coords});
  }
  unless (
    ($pre->{coding}) &&
    $tv->translation_start(undef, $shifting_offset) && 
    $tv->translation_end(undef, $shifting_offset)
  ){
    print "Exiting hgvs_protein - variant " . $vf->variation_name() . "not within translation\n"  if $DEBUG == 1;
    delete($self->{shift_hash}) unless $hash_already_defined;
    return undef;
  }
       
  print "proceeding with hgvs prot\n" if $DEBUG == 1;
  print "Checking translation start: " . $tv->translation_start() ."\n" if $DEBUG == 1;

  ## checks complete - start building term

  ### get reference sequence
  ### this is a temporary fix
  if($tr->stable_id =~ /^ENS|^LRG/ || !defined($tr->analysis) || (defined($tr->analysis) && !defined($tr->analysis->db()))){
    $hgvs_notation->{ref_name} = $tr->translation->display_id();
  }
  else{
    ### get RefSeq identifiers
    my @entries = grep {$_->{dbname} eq 'GenBank'} @{$tr->translation->get_all_DBEntries};
    if(scalar @entries == 1){
      $hgvs_notation->{ref_name} = $entries[0]->{primary_id};
    }
  }

  # Add seq version unless LRG 
  $hgvs_notation->{ref_name} .= "." . $tr->translation->version() 
    unless ($hgvs_notation->{ref_name}=~ /\.\d+$/ || $hgvs_notation->{ref_name} =~ /LRG/);

  $hgvs_notation->{'numbering'} = 'p';

  ### get default reference location [changed later in some cases eg. duplication]
  $hgvs_notation->{start}   = $tv->translation_start();
  $hgvs_notation->{end}     = $tv->translation_end();  

  my $ref = $tv->get_reference_TranscriptVariationAllele;

  ## Incase the user wants shifted HGVS but not shifted consequences, we run the shifting method  
  my $ref_hash_already_defined = defined($ref->{shift_hash});
  $ref->_return_3prime(1) unless $ref_hash_already_defined;
  ## get default reference & alt peptides  [changed later to hgvs format]
  if(defined($self->{shift_hash}) && defined($self->{shift_hash}->{shift_length})  && $self->{shift_hash}->{shift_length} != 0) {
    delete($self->{peptide});
    delete($self->{codon});
    delete($self->{feature_seq});
    $self->shift_feature_seqs();

    if(defined($ref->{shift_hash})) {
      delete($ref->{peptide});
      delete($ref->{codon});
      delete($ref->{feature_seq});

      $ref->{variation_feature_seq} = $self->{shift_hash}->{ref_orig_allele_string};
      $ref->{variation_feature_seq} = $self->{shift_hash}->{shifted_allele_string} if $vf->var_class eq 'deletion';
      $ref->{shifted_feature_seqs} = 1;
    }
  }
  
  $hgvs_notation->{alt} = $self->peptide;

  $hgvs_notation->{ref} = $ref->peptide; 
  
  ## delete the shifting hash if we generated it for HGVS calculations
  delete($self->{shift_hash}) unless $hash_already_defined;
  delete($ref->{shift_hash}) unless $ref_hash_already_defined;

  return undef unless $hgvs_notation->{ref};   
  print "Got protein peps: $hgvs_notation->{ref} =>  $hgvs_notation->{alt} (" . $self->codon() .")\n" if $DEBUG ==1;

  if(defined $hgvs_notation->{alt} && defined $hgvs_notation->{ref} &&
    ($hgvs_notation->{alt} ne  $hgvs_notation->{ref})){
    $hgvs_notation = _clip_alleles( $hgvs_notation);
  }

  #### define type - types are different for protein numbering
  $hgvs_notation  = $self->_get_hgvs_protein_type($hgvs_notation);
  return undef unless defined $hgvs_notation->{type}; 

  ##### Convert ref & alt peptides taking into account HGVS rules
  $hgvs_notation = $self->_get_hgvs_peptides($hgvs_notation);
  unless($hgvs_notation) {
    $self->{hgvs_protein} = undef;
    return undef;
  }

  ##### String formatting
  return $self->_get_hgvs_protein_format($hgvs_notation);
}


=head2 hgvs_offset

  Description: Return the number of bases the variant was shifted 3'
               to defined the HGVS transcript annotation 
  Returntype : int or undef if HGVS has not been calculated or shift not applied
  Exceptions : none
  Status     : At risk

=cut

sub hgvs_offset {
  my $self = shift;
  #_hgvs_offset can usually be taken directly from the shift hash, however in situations where we remove the shift_hash from the $tva object after calculating HGVS then we can access it from $self->{_hgvs_offset}
  return defined($self->{shift_hash}) ? $self->{shift_hash}->{_hgvs_offset} : $self->{_hgvs_offset};
}

=head2 hgvs_exon_start_coordinate

  Description: Return the HGVS exon start coordinate
  Returntype : int or undef if HGVS has not been calculated
  Exceptions : none
  Status     : At risk

=cut

sub hgvs_exon_start_coordinate {
  my $self = shift;
  return $self->{_hgvs_exon_start_coordinate};
}

=head2 hgvs_intron_start_offset

  Description: Return the HGVS intron start offset
  Returntype : int or undef if HGVS has not been calculated
  Exceptions : none
  Status     : At risk

=cut

sub hgvs_intron_start_offset {
  my $self = shift;
  return $self->{_hgvs_intron_start_offset};
}

=head2 hgvs_exon_end_coordinate

  Description: Return the HGVS exon end coordinate
  Returntype : int or undef if HGVS has not been calculated
  Exceptions : none
  Status     : At risk

=cut

sub hgvs_exon_end_coordinate {
  my $self = shift;
  return $self->{_hgvs_exon_end_coordinate};
}

=head2 hgvs_intron_end_offset

  Description: Return the HGVS intron end offset
  Returntype : int or undef if HGVS has not been calculated
  Exceptions : none
  Status     : At risk

=cut

sub hgvs_intron_end_offset {
  my $self = shift;
  return $self->{_hgvs_intron_end_offset};
}

### HGVS: format protein string
sub _get_hgvs_protein_format {
  my $self          = shift;
  my $hgvs_notation = shift;

  ### all start with refseq name & numbering type
  $hgvs_notation->{'hgvs'} = $hgvs_notation->{'ref_name'} . ":" . $hgvs_notation->{'numbering'} . ".";    

  ### New (v 15.11) way to describe synonymous changes
  if( $hgvs_notation->{ref} eq $hgvs_notation->{alt} 
       && $hgvs_notation->{type} ne "fs" && $hgvs_notation->{type} ne "ins"){
    return $hgvs_notation->{'hgvs'} . $hgvs_notation->{ref} . $hgvs_notation->{start} . "=";
  }

  ### handle stop_lost seperately regardless of cause by del/delins => p.TerposAA1extnum_AA_to_stop
  if(stop_lost($self) && ($hgvs_notation->{type} eq "del" || $hgvs_notation->{type} eq ">" )) {
    ### if deletion of stop add extTer and number of new aa to alt

    $hgvs_notation->{alt} = substr($hgvs_notation->{alt}, 0, 3);
    print "stop loss check req for $hgvs_notation->{type}\n" if $DEBUG ==1;

    my $aa_til_stop =  $self->_stop_loss_extra_AA($hgvs_notation->{start}-1 );
    ### use ? to show new stop not predicted
    $aa_til_stop = "?" unless defined $aa_til_stop ;

    $hgvs_notation->{alt} .=  "extTer" . $aa_til_stop;

    # 2020-07-14
    if( length($hgvs_notation->{ref}) >3 && $hgvs_notation->{type} eq "del" ) {
      my $ref_pep_first = substr($hgvs_notation->{ref}, 0, 3);
      my $ref_pep_last  = substr($hgvs_notation->{ref}, -3, 3);
      $hgvs_notation->{'hgvs'} .=  $ref_pep_first . $hgvs_notation->{start} .  "_" .  $ref_pep_last . $hgvs_notation->{end} .$hgvs_notation->{alt} ;
    }
    else{
      $hgvs_notation->{'hgvs'} .=  $hgvs_notation->{ref} . $hgvs_notation->{start} .  $hgvs_notation->{alt} ;
    }

  } 

  elsif( $hgvs_notation->{type} eq "dup"){

    if($hgvs_notation->{start} < $hgvs_notation->{end}){
      ### list only first and last peptides in long duplicated string
      my $ref_pep_first = substr($hgvs_notation->{alt}, 0, 3);
      my $ref_pep_last  = substr($hgvs_notation->{alt}, -3, 3);
      $hgvs_notation->{'hgvs'} .=  $ref_pep_first . $hgvs_notation->{start} .  "_" .  $ref_pep_last . $hgvs_notation->{end} ."dup";
    }
    else{
      $hgvs_notation->{'hgvs'} .=  $hgvs_notation->{alt} . $hgvs_notation->{start} .  "dup" ;
    }

    print "formating a dup  $hgvs_notation->{hgvs} \n" if $DEBUG==1;
  }

  elsif($hgvs_notation->{type} eq ">"){
    #### substitution
    $hgvs_notation->{'hgvs'}  .=   $hgvs_notation->{ref}. $hgvs_notation->{start} .  $hgvs_notation->{alt};
  }    

  elsif( $hgvs_notation->{type} eq "delins" || $hgvs_notation->{type} eq "ins" ){

    ## don't report other peptides after a stop is found
    $hgvs_notation->{alt} =~ s/Ter\w+/Ter/ ;

    #### list first and last aa in reference only
    my $ref_pep_first = substr($hgvs_notation->{ref}, 0, 3);
    my $ref_pep_last;
    if(substr($hgvs_notation->{ref}, -1, 1) eq "X"){
      $ref_pep_last ="Ter";
    }
    else{
      $ref_pep_last  = substr($hgvs_notation->{ref}, -3, 3);
    }

    if($hgvs_notation->{ref} =~ /X$/) {
      ### For stops & add extX & distance to next stop to alt pep
      my $aa_til_stop =  $self->_stop_loss_extra_AA( $hgvs_notation->{start}-1, "loss");
      if(defined  $aa_til_stop){  
         $hgvs_notation->{alt} .="extTer" . $aa_til_stop;
      }
    }


    if($hgvs_notation->{start} == $hgvs_notation->{end} && $hgvs_notation->{type} eq "delins"){       
       $hgvs_notation->{'hgvs'} .= $ref_pep_first . $hgvs_notation->{start}  . $hgvs_notation->{type} . $hgvs_notation->{alt} ;
    }
    else{        
      ### correct ordering if needed
      if($hgvs_notation->{start} > $hgvs_notation->{end}){       
        ( $hgvs_notation->{start}, $hgvs_notation->{end}) = ($hgvs_notation->{end}, $hgvs_notation->{start} );
      }

      $hgvs_notation->{'hgvs'} .= $ref_pep_first . $hgvs_notation->{start}  . "_"  .  $ref_pep_last . $hgvs_notation->{end} . $hgvs_notation->{type} . $hgvs_notation->{alt} ;
    }
  }     
  
  elsif($hgvs_notation->{type} eq "fs"){

    if(defined $hgvs_notation->{alt} && $hgvs_notation->{alt} eq "Ter"){ ## stop gained
      ## describe as substitution if stop occurs immediately
      $hgvs_notation->{'hgvs'} .= $hgvs_notation->{ref} . $hgvs_notation->{start}  .  $hgvs_notation->{alt} ;
    }
    else{ 
      ## not immediate stop - count aa until next stop
      my $aa_til_stop =  $self->_stop_loss_extra_AA( $hgvs_notation->{start}-1, "fs");

      ### use ? to show new stop not predicted
      $aa_til_stop = "?" unless defined $aa_til_stop; 

      $hgvs_notation->{'hgvs'} .= $hgvs_notation->{ref} . $hgvs_notation->{start} . $hgvs_notation->{alt} ."fsTer$aa_til_stop";     
    }
  }

  elsif( $hgvs_notation->{type} eq "del"){
    $hgvs_notation->{alt} =  "del"; 
    if( length($hgvs_notation->{ref}) >3 ){
      my $ref_pep_first = substr($hgvs_notation->{ref}, 0, 3);
      my $ref_pep_last  = substr($hgvs_notation->{ref}, -3, 3);
      $hgvs_notation->{'hgvs'} .=  $ref_pep_first . $hgvs_notation->{start} .  "_" .  $ref_pep_last . $hgvs_notation->{end} .$hgvs_notation->{alt} ;
    }
    else{
      $hgvs_notation->{'hgvs'} .=  $hgvs_notation->{ref} . $hgvs_notation->{start} .  $hgvs_notation->{alt} ;
    }       
  }

  elsif($hgvs_notation->{start} ne $hgvs_notation->{end} ){
    $hgvs_notation->{'hgvs'}  .=  $hgvs_notation->{ref} . $hgvs_notation->{start}  . "_" .  $hgvs_notation->{alt} . $hgvs_notation->{end} ;
  }

  else{
  #### default to substitution    
    $hgvs_notation->{'hgvs'}  .=   $hgvs_notation->{ref}. $hgvs_notation->{start} .  $hgvs_notation->{alt};
  }

  if($DEBUG==1){ print "Returning protein format: $hgvs_notation->{'hgvs'}\n";}
  return $hgvs_notation->{'hgvs'};
}

### HGVS: get type of variation event in protein terms
sub _get_hgvs_protein_type {
  my $self = shift;
  my $hgvs_notation = shift;

  if($DEBUG==1){ print "starting get_hgvs_protein_type \n";}

  if( frameshift($self) ){
    $hgvs_notation->{type} = "fs";
    return $hgvs_notation;
  }

  if( defined $hgvs_notation->{ref} && defined $hgvs_notation->{alt} ){
  ### Run type checks on peptides if available
    $hgvs_notation->{ref} =~ s/\*/X/;
    $hgvs_notation->{alt} =~ s/\*/X/;

    if($hgvs_notation->{ref} eq "-" || $hgvs_notation->{ref} eq "") {
      $hgvs_notation->{type} = "ins";
    }
    elsif($hgvs_notation->{alt} eq "" || $hgvs_notation->{alt} eq "-") {
      $hgvs_notation->{type} = "del";
    }
    elsif( length($hgvs_notation->{ref}) ==1 && length($hgvs_notation->{alt}) ==1 ) {
      $hgvs_notation->{type} = ">";
    }
    elsif(
      ((length($hgvs_notation->{alt}) >0 && length($hgvs_notation->{ref}) >0) &&
      (length($hgvs_notation->{alt}) ne length($hgvs_notation->{ref})) )  ||
      (length($hgvs_notation->{alt}) >1 && length($hgvs_notation->{ref}) >1)     ## not a substitution if >1 aa switched
    ) {
      $hgvs_notation->{type} = "delins";
    }
    else{
      $hgvs_notation->{type} = ">";
    }
  }
  else {
    ### Cannot define type from peptides - check at DNA level
    ### get allele length from dna seq & cds length
    my ($ref_length, $alt_length ) = $self->_get_allele_length(); 

    if($alt_length >1  ){
      if($hgvs_notation->{start} == ($hgvs_notation->{end} + 1) ){
        ### convention for insertions - end one less than start
        $hgvs_notation->{type} = "ins";
      }
      elsif( $hgvs_notation->{start} != $hgvs_notation->{end}  ){
        $hgvs_notation->{type} = "delins";
      }
      else{
        $hgvs_notation->{type} = ">";
      } 
    }
    elsif($ref_length >1  ){
      $hgvs_notation->{type} = "del"; 
    }
 
    else{
      #print STDERR "DEBUG ".$self->variation_feature->start."\n";
      #warn "Cannot define protein variant type [$ref_length  - $alt_length]\n";
    }
  }

  return $hgvs_notation;
}

### HGVS: get reference & alternative peptide 
sub _get_hgvs_peptides {
  my $self          = shift;
  my $hgvs_notation = shift;

  if($hgvs_notation->{type} eq "fs"){
    ### ensembl alt/ref peptides not the same as HGVS alt/ref - look up seperately
    $hgvs_notation = $self->_get_fs_peptides($hgvs_notation); 
    return undef unless defined $hgvs_notation->{type};   
  }
  elsif($hgvs_notation->{type} eq "ins" ){

    ### Check if bases directly after insertion match inserted sequence
    $hgvs_notation = $self->_check_peptides_post_var($hgvs_notation);

    ### Check that inserted bases do not duplicate 3' reference sequence [set to type = dup and return if so]
    $hgvs_notation = $self->_check_for_peptide_duplication($hgvs_notation) unless $hgvs_notation->{alt} =~/\*/;
    return ($hgvs_notation) if $hgvs_notation->{type} eq "dup";              

    ### HGVS ref are peptides flanking insertion
    my $min;
    if($hgvs_notation->{start} < $hgvs_notation->{end}){
        $min = $hgvs_notation->{start};
    }
    else{ $min = $hgvs_notation->{end};}

    $hgvs_notation->{ref} = $self->_get_surrounding_peptides(
      $min, 
      $hgvs_notation->{original_ref},
      2
    );

    return undef unless $hgvs_notation->{ref};
  }
  elsif($hgvs_notation->{type} eq "del" ){
    ##check if bases directly after deletion match the deletion
    $hgvs_notation = $self->_check_peptides_post_var($hgvs_notation);
  }

  ### Convert peptide to 3 letter code as used in HGVS
  unless( $hgvs_notation->{ref} eq "-"){
    $hgvs_notation->{ref}  = Bio::SeqUtils->seq3(Bio::PrimarySeq->new(-seq => $hgvs_notation->{ref}, -id => 'ref',  -alphabet => 'protein')) || "";
  }
  if( $hgvs_notation->{alt} eq "-"){
    $hgvs_notation->{alt} = "del";
  }
  else{
    $hgvs_notation->{alt}  = Bio::SeqUtils->seq3(Bio::PrimarySeq->new(-seq => $hgvs_notation->{alt}, -id => 'ref',  -alphabet => 'protein')) || "";
  }

  ### handle special cases
  if( start_lost($self) ){
    #### handle initiator loss -> probably no translation => alt allele is '?'
    $hgvs_notation->{alt}  = "?";    
    $hgvs_notation->{type} = "";
  }

  elsif(  $hgvs_notation->{type} eq "del"){ 
    if( $hgvs_notation->{ref} =~/\w+/){
      $hgvs_notation->{alt} = "del";
    }
    else{
      $hgvs_notation = $self->_get_del_peptides($hgvs_notation);
    }
  } 
  elsif($hgvs_notation->{type} eq "fs"){
    ### only quote first ref peptide for frameshift
    $hgvs_notation->{ref} = substr($hgvs_notation->{ref},0,3);
  }

  ### 2012-08-31  - Ter now recommended in HGVS
  if(defined $hgvs_notation->{ref}){ $hgvs_notation->{ref} =~ s/Xaa/Ter/g; }
  if(defined $hgvs_notation->{alt}){ $hgvs_notation->{alt} =~ s/Xaa/Ter/g; }

  return ($hgvs_notation);           
}

### HGVS:  remove common peptides/nucleotides from alt and ref strings & shift start/end accordingly 
sub _clip_alleles {
    
  my $hgvs_notation = shift;
    
  my $check_alt   = $hgvs_notation->{alt};
  my $check_ref   = $hgvs_notation->{ref};
  my $check_start = $hgvs_notation->{start};
  my $check_end   = $hgvs_notation->{end};

  ## cache this - if stop needed later
  $hgvs_notation->{original_ref} = $hgvs_notation->{ref};

  ## store identical trimmed seq 
  my $preseq;
  print "can we clip :  $check_ref &  $check_alt\n" if $DEBUG ==1;
  ### strip same bases from start of string
  for (my $p =0; $p <length ($hgvs_notation->{ref}); $p++){
    my $check_next_ref = substr( $check_ref, 0, 1);
    my $check_next_alt = substr( $check_alt, 0, 1);
    
    if( defined $hgvs_notation->{'numbering'} && $hgvs_notation->{'numbering'} eq 'p' && 
        $check_next_ref eq  "*" && $check_next_alt eq "*"){
      ### stop re-created by variant - no protein change
      $hgvs_notation->{type} = "=";

      return($hgvs_notation);
    }
    
    if($check_next_ref eq  $check_next_alt){
      $check_start++;
      $check_ref  = substr( $check_ref, 1);
      $check_alt  = substr( $check_alt, 1);
      $preseq    .= $check_next_ref;
    }
    else{
      last;
    }
  }

  my $len = length ($check_ref);
  #### strip same bases from end of string

  for (my $q =0; $q < $len; $q++) {
    my $check_next_ref = substr( $check_ref, -1, 1);
    my $check_next_alt = substr( $check_alt, -1, 1);
    if($check_next_ref eq  $check_next_alt){
      chop $check_ref;
      chop $check_alt;
      $check_end--;
    }
    else{
      last;
    }
  }
  
  ## ammend positions & ref/alt
  $hgvs_notation->{alt}   = $check_alt;
  $hgvs_notation->{ref}   = $check_ref;
  $hgvs_notation->{start} = $check_start;
  $hgvs_notation->{end}   = $check_end;
  $hgvs_notation->{preseq} =   $preseq ;

  ### check if clipping suggests a type change 

  ## no protein change - use transcript level annotation 
  $hgvs_notation->{type} = "="   if( defined $hgvs_notation->{'numbering'} && 
                                     $hgvs_notation->{'numbering'} eq 'p' &&
                                     $hgvs_notation->{alt} eq $hgvs_notation->{ref});      

  ### re-set as ins not delins    
  $hgvs_notation->{type} ="ins"  if(length ($check_ref) == 0 && length ($check_alt) >= 1);

  ### re-set as del not delins  
  $hgvs_notation->{type}  ="del" if(length ($check_ref) >=1 && length ($check_alt) == 0);      

  print "clipped :  $check_ref &  $check_alt\n" if $DEBUG ==1;

  return $hgvs_notation;
}
    




#### HGVS: check allele lengths to look for frameshifts
sub _get_allele_length {
  my $self       = shift;
  my $ref_length = 0;
  my $alt_length = 0;

  my $al_string = $self->allele_string();
  my $ref_allele = (split/\//, $al_string)[0];
  $ref_allele =~ s/\-//;
  $ref_length = length $ref_allele;

  my $alt_allele = $self->variation_feature_seq();
  $alt_allele =~ s/\-//;
  $alt_length = length $alt_allele;

  return ($ref_length, $alt_length );  
}

### HGVS: list first different peptide [may not be first changed codon]
sub _get_fs_peptides {
  my $self    = shift;
  my $hgvs_notation = shift;

  ### get CDS with alt variant
  my $alt_cds = $self->_get_alternate_cds();
  return undef unless defined($alt_cds);
  unless ($alt_cds->seq =~ /([ACGT-])+/) {
    warn('Alternate CDS not found on seq region ' . $self->variation_feature->seq_region_name . '. Are you missing a synonyms file?');
    return undef;
  }

  #### get new translation
  my $alt_trans = $alt_cds->translate()->seq();

  ### get changed end (currently in single letter AA coding)    
  my $ref_trans  = $self->transcript_variation->_peptide;

  $ref_trans    .= "*";   ## appending ref stop for checking purposes 
  
  $hgvs_notation->{start} = $self->transcript_variation->translation_start() ;

  if( $hgvs_notation->{start} > length $alt_trans){ ## deletion of stop, no further AA in alt seq 
    $hgvs_notation->{alt}  = "del";
    $hgvs_notation->{type} = "del";
    return $hgvs_notation;
  }

  while ($hgvs_notation->{start} <= length $alt_trans){
    ### frame shift may result in the same AA#

    $hgvs_notation->{ref} = substr($ref_trans, $hgvs_notation->{start}-1, 1);
    $hgvs_notation->{alt} = substr($alt_trans, $hgvs_notation->{start}-1, 1);

    if($hgvs_notation->{ref} eq "*" && $hgvs_notation->{alt} eq "*"){
      ### variation at stop codon, but maintains stop codon - set to synonymous
      $hgvs_notation->{type} = "=";
      return ($hgvs_notation);
    }

    last if $hgvs_notation->{ref} ne $hgvs_notation->{alt};
    $hgvs_notation->{start}++;
  }

  return ($hgvs_notation);
}

#### HGVS: if variant is an insertion, ref pep is initially "-", so seek AA before and after insertion
sub _get_surrounding_peptides {
  my $self    = shift;
  my $ref_pos = shift; 
  my $original_ref = shift;
  my $length  = shift;

  my $ref_trans  = $self->transcript_variation->_peptide();
  $ref_trans .= $original_ref
    if defined $original_ref && $original_ref =~ /^\*/;

  ## can't find peptide after the end
  return if length($ref_trans) <=  $ref_pos ;

  my $ref_string;
  if(defined $length) {
    $ref_string = substr($ref_trans, $ref_pos-1, $length);
  }
  else{
    $ref_string = substr($ref_trans, $ref_pos-1 );
  }

  return ($ref_string);
}


#### HGVS: alternate CDS needed to check for new stop when variant disrupts 'reference' stop
sub _get_alternate_cds {
    
  my $self = shift;

  ### get reference sequence
  my $reference_cds_seq = $self->transcript_variation->_translateable_seq();
  
  my $tv = $self->transcript_variation;
  my $vf = $tv->variation_feature;
  my $tr = $tv->transcript;
  my $shifting_offset = (defined($self->{shift_hash}) && defined($self->{shift_hash}->{shift_length})) ? $self->{shift_hash}->{shift_length} : 0;
  
  return undef unless defined($tv->cds_start(undef, $tr->strand * $shifting_offset)) && defined($tv->cds_end(undef, $tr->strand * $shifting_offset));

  ### get sequences upstream and downstream of variant
  my $upstream_seq   =  substr($reference_cds_seq, 0, ($tv->cds_start(undef, $tr->strand * $shifting_offset) -1) );
  my $downstream_seq =  substr($reference_cds_seq, ($tv->cds_end(undef, $tr->strand * $shifting_offset)) );
  return undef unless defined($downstream_seq) && defined($upstream_seq);
  
  ### fix alternate allele if deletion or on opposite strand
  my $alt_allele  = $self->variation_feature_seq();
  $alt_allele  =~ s/\-//;
  if($alt_allele && $vf->strand() != $tr->strand()){    
    reverse_comp(\$alt_allele) ;
  }

  ### build alternate seq
  my $alternate_seq  = $upstream_seq . $alt_allele . $downstream_seq ;
  $alternate_seq  = $self->_trim_incomplete_codon($alternate_seq );

  ### create seq obj with alternative allele in the CDS sequence
  my $alt_cds =Bio::PrimarySeq->new(-seq => $alternate_seq,  -id => 'alt_cds', -alphabet => 'dna');

  ### append UTR if available as stop may be disrupted
  my $utr = $self->transcript_variation->_three_prime_utr();

  if (defined $utr) {
  ### append the UTR to the alternative CDS 
    $alt_cds->seq($alt_cds->seq() . $utr->seq()); 
  }
  else{
   ##warn "No UTR available for alternate CDS\n";
  }

  return $alt_cds;
}

### HGVS: if inserted string is identical to 3' reference sequence, describe as duplication
sub _check_for_peptide_duplication {    
  my $self = shift;
  my $hgvs_notation = shift;

  ##### get reference sequence
  my $reference_cds_seq = $self->transcript_variation->_translateable_seq();

  my $reference_cds = Bio::PrimarySeq->new(-seq => $reference_cds_seq,  -id => 'alt_cds', -alphabet => 'dna');
  my $reference_trans = $reference_cds->translate()->seq();

  ##### get sequence upstream of variant - use hgvs start; may have been shifted
  my $upstream_trans  = substr($reference_trans, 0, ($hgvs_notation->{'start'} -1) );
  print "Checking for peptide duplication: $hgvs_notation->{alt} vs $upstream_trans  $hgvs_notation->{preseq} \n" if $DEBUG ==1;

  $upstream_trans  .= $hgvs_notation->{preseq} if defined $hgvs_notation->{preseq}; ## add back on anything previously chopped off ref allele

  ## Test whether alt peptide matches the reference sequence just before the variant
  my $test_new_start = $hgvs_notation->{'start'} - length($hgvs_notation->{'alt'}) -1 ;

  if( (length($upstream_trans) >=  $test_new_start + length($hgvs_notation->{'alt'}) ) && $test_new_start  >=0){
    my $test_seq       =  substr($upstream_trans, $test_new_start, length($hgvs_notation->{'alt'}));

    if ( $test_new_start >= 0 && $test_seq eq $hgvs_notation->{alt}) {           
      $hgvs_notation->{type}   = 'dup';
      $hgvs_notation->{end}    = $hgvs_notation->{start} -1;
      $hgvs_notation->{start} -= length($hgvs_notation->{alt});

      ## convert to 3 letter code
      $hgvs_notation->{alt}  = Bio::SeqUtils->seq3(Bio::PrimarySeq->new(-seq => $hgvs_notation->{alt}, -id => 'ref',  -alphabet => 'protein')) || "";
    }
  }
  return $hgvs_notation;
}

#### HGVS: if a stop is lost, seek the next in transcript & count number of extra AA's
sub _stop_loss_extra_AA{

  my $self        = shift;
  my $ref_var_pos = shift;  ### first effected AA - supply for frameshifts
  my $test        = shift;

  return undef unless $ref_var_pos;

  my $extra_aa;

  ### get the sequence with variant added
  my $alt_cds   = $self->_get_alternate_cds();
  return undef unless defined($alt_cds);
  
  ### get new translation
  my $alt_trans = $alt_cds->translate();

  my $ref_temp  =  $self->transcript_variation->_peptide();
  my $ref_len = length($ref_temp);
  
  if($DEBUG==1){ 
    print "alt translated:\n" . $alt_trans->seq() . "\n";
    print "ref translated:\n$ref_temp\n";;
  }
  
  #### Find the number of residues that are translated until a termination codon is encountered
  if ($alt_trans->seq() =~ m/\*/) {
    if($DEBUG==1){print "Got $+[0] aa before stop, var event at $ref_var_pos \n";}
  
    if(defined $test && $test eq "fs" ){
      ### frame shift - count from first AA effected by variant to stop
      $extra_aa = $+[0] - $ref_var_pos;
      if($DEBUG==1){ print "Stop change ($test): found $extra_aa amino acids before fs stop [ $+[0] - peptide ref_start: $ref_var_pos )]\n";}
    }
  
    else{
      $extra_aa = $+[0]  - 1 - $ref_len;
      if($DEBUG==1){ print "Stop change (non-fs): found $extra_aa amino acids before next stop [ $+[0] - 1 -normal stop $ref_len)]\n";}        
    }
  }
  
  # A special case is if the first aa is a stop codon => don't display the number of residues until the stop codon
  if(defined $extra_aa && $extra_aa >0){ 
    return $extra_aa ;
  }
  else{ 
    #warn "No stop found in alternate sequence\n";
    return undef;
  }
}

## doing this to stop final incomplete codons being guessed
sub _trim_incomplete_codon{
  my $self = shift;
  my $seq  = shift;

  return undef unless $seq;

  my $full_length = length $seq;
  my $keep_length = $full_length -  $full_length % 3;
  return $seq if $full_length = $keep_length ;

  return substr($seq, 0, $keep_length);
}


## This is used for rare in-frame deletions removing an intron and part of both surrounding exons
sub _get_del_peptides{

  my $self    = shift;
  my $hgvs_notation = shift;

  ### get CDS with alt variant
  my $alt_cds = $self->_get_alternate_cds();
  return undef unless defined($alt_cds);

  #### get new translation
  my $start = $self->transcript_variation->translation_start() - 1;
  my $alt = substr($alt_cds->translate()->seq(), $start );
  $hgvs_notation->{alt} = (split/\*/, $alt)[0];

  ### get changed end (currently in single letter AA coding)    
  $hgvs_notation->{ref}  = substr($self->transcript->translate()->seq(), $start );

  $hgvs_notation->{start} = $self->transcript_variation->translation_start() ;

  $hgvs_notation = _clip_alleles($hgvs_notation);

  ## switch to 3 letter AA coding
  $hgvs_notation->{alt}  = Bio::SeqUtils->seq3(Bio::PrimarySeq->new(-seq => $hgvs_notation->{alt}, -id => 'ref',  -alphabet => 'protein')) || "";
  $hgvs_notation->{ref}  = Bio::SeqUtils->seq3(Bio::PrimarySeq->new(-seq => $hgvs_notation->{ref}, -id => 'ref',  -alphabet => 'protein')) || "";
  return $hgvs_notation;
}

## HGVS counts from first different peptide, 
## so check the sequence post variant and increment accordingly
sub _check_peptides_post_var{
  my $self          = shift;
  my $hgvs_notation = shift;

  ## check peptides after deletion 
  my $post_pos = $hgvs_notation->{end}+1;
  my $post_seq = $self->_get_surrounding_peptides(
    $post_pos,
    $hgvs_notation->{original_ref}
  );

  ## if a stop is deleted and no sequence is available beyond to check, return
  return $hgvs_notation unless defined $post_seq;

  $hgvs_notation = _shift_3prime($hgvs_notation, $post_seq);
  return $hgvs_notation;
}

## HGVS aligns changes 3' (alt string may look different at transcript level to genomic level)
## AAC[TG]TAT => AACT[GT]AT
## TTA[GGG]GGTTTA =>TTAGG[GGG]TTTA

sub _shift_3prime{

  my $hgvs_notation = shift;
  my $post_seq      = shift;

  my $seq_to_check;
  if( $hgvs_notation->{type} eq 'ins'){
    $seq_to_check = $hgvs_notation->{alt};
  }
  elsif ($hgvs_notation->{type} eq 'del'){
    $seq_to_check = $hgvs_notation->{ref};
  }
  else{
    return $hgvs_notation;
  }

  ## return if nothing to check
  return $hgvs_notation unless defined $post_seq && defined $seq_to_check;

  ## get length of pattern to check 
  my $deleted_length = (length $seq_to_check);
  
  # warn "Checking $seq_to_check v $post_seq\n";
  ## move along sequence after deletion looking for match to start of deletion
  for (my $n = 0; $n<= (length($post_seq) - $deleted_length); $n++ ){

    ## check each position in deletion/ following seq for match
    my $check_next_del  = substr( $seq_to_check, 0, 1);
    my $check_next_post = substr( $post_seq, $n, 1);

    if($check_next_del eq $check_next_post){

      ## move position of deletion along
      $hgvs_notation->{start}++;
      $hgvs_notation->{end}++;
  	
      ## modify deleted sequence - remove start & append to end
      $seq_to_check = substr($seq_to_check,1);
      $seq_to_check .= $check_next_del;
    }
    else{
      last;	    
    }
  }
  ## set new HGVS string
  $hgvs_notation->{alt} = $seq_to_check if $hgvs_notation->{type} eq 'ins';
  $hgvs_notation->{ref} = $seq_to_check if $hgvs_notation->{type} eq 'del';

  return $hgvs_notation;
}
=head
# We haven't implemented support for these methods yet

sub hgvs_rna {
    return _hgvs_generic(@_,'rna');
}

sub hgvs_mitochondrial {
    return _hgvs_generic(@_,'mitochondrial');
}

=cut

sub _hgvs_generic {
  my $self = shift;
  my $reference = pop;
  my $notation = shift;
  
  #The rna and mitochondrial modes have not yet been implemented, so return undef in case we get a call to these
  return undef if ($reference =~ m/rna|mitochondrial/);
  
  my $sub = qq{hgvs_$reference};
  
  $self->{$sub} = $notation if defined $notation;
  
  unless ($self->{$sub}) {

    # Use the transcript this VF is on as the reference feature
    my $reference_feature = $self->feature;
    # If we want genomic coordinates, the reference_feature should actually be the slice for the underlying seq_region
    $reference_feature = $reference_feature->slice->seq_region_Slice if ($reference eq 'genomic');

    # Calculate the HGVS notation on-the-fly and pass it to the TranscriptVariation in order to distribute the result to the other alleles
    my $tv = $self->base_variation_feature_overlap;
    my $vf = $self->base_variation_feature;

    $tv->$sub($vf->get_all_hgvs_notations($reference_feature,substr($reference,0,1),undef,undef,$tv));
  }
  
  return $self->{$sub};
}


### HGVS: move variant to transcript slice
sub _var2transcript_slice_coords{
  my ($self, $tr, $tv, $vf) = @_;

  $tv ||= $self->base_variation_feature_overlap;
  $tr ||= $tv->transcript;
  $vf ||= $tv->base_variation_feature;
  
  # what we want is VF coords relative to transcript feature slice
  my ($tr_start, $tr_end) = ($tr->start, $tr->end);
  my ($vf_start, $vf_end);

  # same slice, easy and we don't need to use transfer
  if($NO_TRANSFER || $tr->slice eq $vf->slice) {

    # different transform depending on transcript strand
    if($tr->strand < 1) {

      # note we switch start/end here
      # this also works in the case of insertions thankfully
      ($vf_start, $vf_end) = map {($tr_end - $_) + 1} ($vf->end, $vf->start) unless $vf->{shifted_flag};
      ($vf_start, $vf_end) = map {($tr_end - $_) + 1} ($vf->{unshifted_end}, $vf->{unshifted_start}) if $vf->{shifted_flag};
    }
    else {
      ($vf_start, $vf_end) = map {($_ - $tr_start) + 1} ($vf->start, $vf->end) unless $vf->{shifted_flag};
      ($vf_start, $vf_end) = map {($_ - $tr_start) + 1} ($vf->{unshifted_start}, $vf->{unshifted_end}) if $vf->{shifted_flag};
    }
  }

  # different slices used to fetch features
  # have to use transfer for safety
  else {
    my $tr_vf = $vf->transfer($self->_transcript_feature_Slice($tr));
    return undef unless $tr_vf;
    ($vf_start, $vf_end) = ($tr_vf->start, $tr_vf->end);
  }

  # Return undef if this VariationFeature does not fall within the supplied feature.
  return undef if (
    $vf_start  < 1 || 
    $vf_end    < 1 || 
    $vf_start  > ($tr_end - $tr_start + 1) || 
    $vf_end    > ($tr_end - $tr_start + 1)
  ); 
  
  return( $vf_start , $vf_end, $self->_transcript_feature_Slice($tr));
}



### HGVS: get variant position in transcript 

# intron: 
# If the position is in an intron, the boundary position of the closest exon and 
# a + or - offset into the intron is returned.
# Ordered by genome forward not 5' -> 3'

# upstream:
# If the position is 5' of the start codon, it is reported relative to the start codon 
# (-1 being the last nucleotide before the 'A' of ATG).

#downstream:
# If the position is 3' pf the stop codon, it is reported with a '*' prefix and the offset 
# from the start codon (*1 being the first nucleotide after the last position of the stop codon)

sub _get_cDNA_position {

  my $self     = shift;
  my $position = shift; ### start or end of variation

  my $tv         = $self->base_variation_feature_overlap;
  my $transcript = $tv->transcript();
  my $strand     = $transcript->strand();    

  #### TranscriptVariation start/stop coord relative to transcript 
  #### Switch to chromosome coordinates taking into account strand
  $position = ( $strand > 0 ? 
    ( $transcript->start() + $position - 1 )  :   
    ( $transcript->end()   - $position + 1));

  # Get all exons sorted in positional order
  my $exons = $tv->_sorted_exons();

  my $n_exons = scalar(@$exons);

  my $cdna_position;
  # Loop over the exons and get the coordinates of the variation in exon+intron notation
  for (my $i=0; $i<$n_exons; $i++) {

    my $exon = $exons->[$i];
    my ($exon_start, $exon_end) = ($exon->{start}, $exon->{end});

    # Skip if the start point is beyond this exon
    next if ($position > $exon_end);


    # EXONIC:  If the start coordinate is within this exon
    if ($position >= $exon_start) {
      # Get the cDNA start coordinate of the exon and add the number of nucleotides from the exon boundary to the variation
      # If the transcript is in the opposite direction, count from the end instead
      $cdna_position = $self->_exon_cdna_start($exon, $transcript) + (
        $strand > 0 ? 
        ( $position - $exon_start ) : 
        ( $exon_end - $position ) 
      );
      last;  #### last exon checked
    }

    ## INTRONIC
    # Else the start coordinate is between this exon and the previous one, determine which one is closest and get coordinates relative to that one
    else {

      my $prev_exon = $exons->[$i-1];

      my $updist   = ($position - $prev_exon->{end});
      $updist =~ s/\-//; ## avoid problems with incomplete transcripts
      my $downdist = ($exon_start - $position);
      $downdist =~ s/\-//; ## avoid problems with incomplete transcripts

      # If the distance to the upstream exon is the shortest, or equal and in the positive orientation, use that
      if ($updist < $downdist || ($updist == $downdist && $strand >= 0)) {
        
        # If the orientation is reversed, we should use the cDNA start and a '-' offset
        $cdna_position = (
          $strand >= 0 ? 
          $self->_exon_cdna_end($prev_exon, $transcript).'+' : 
          $self->_exon_cdna_start($prev_exon, $transcript).'-'
        ).$updist;
      }
      # Else if downstream is shortest...
      else {
        # If the orientation is reversed, we should use the cDNA end and a '+' offset
        $cdna_position = (
          $strand >= 0 ?
          $self->_exon_cdna_start($exon, $transcript).'-' : 
          $self->_exon_cdna_end($exon, $transcript).'+'
        ).$downdist;
      }

      last; ## last exon checked
    }
  }

  ## this should not happen; potential refseq oddness
  return undef unless $cdna_position; 
 
  # Shift the position to make it relative to the start & stop codons
  my $start_codon  =  $transcript->cdna_coding_start();
  my $stop_codon   =  $transcript->cdna_coding_end();

  # Disassemble the cDNA coordinate into the exon and intron parts
  ### just built this now taking it appart again
  my ($cdna_coord, $intron_offset) = $cdna_position =~ m/([0-9]+)([\+\-][0-9]+)?/;
  

  # Start by correcting for the stop codon
  if (defined($stop_codon) ){

    if($cdna_coord > $stop_codon) {
      # Get the offset from the stop codon
      $cdna_coord -= $stop_codon;
      # Prepend a * to indicate the position is in the 3' UTR
      $cdna_coord = '*' . $cdna_coord;
    }
    elsif ( $cdna_coord eq $stop_codon && defined $intron_offset) {
      $intron_offset =~ s/\+//g;
      $cdna_coord ='' ;

      # Prepend a * to indicate the position is in the 3' UTR
      $cdna_coord = '*' . $cdna_coord;
    }
  }
  if (defined($start_codon) && $cdna_coord  !~/\*/) {
    
    # If the position is beyond the start codon, add 1 to get the correct offset
    $cdna_coord += ($cdna_coord >= $start_codon);
    # Subtract the position of the start codon
    $cdna_coord -= $start_codon;
  }
  else{
  	print "Checking non-coding transcript\n" if $DEBUG==1;
  }

  # Re-assemble the cDNA position  [ return exon num & offset & direction for intron eg. 142+363] 
  $cdna_position = $cdna_coord . (defined($intron_offset) ? $intron_offset : '');

  return $cdna_position;
}

# $exon->cdna_start doesn't cache
# so use our own method that does
sub _exon_cdna_start {
  my ($self, $exon, $tr) = @_;

  my $tr_stable_id = $tr->stable_id;
  my $fc = $exon->{_variation_effect_feature_cache}->{$tr_stable_id} ||= {};

  if(!exists($fc->{_cdna_start})) {
    $fc->{_cdna_start} = $exon->cdna_start($tr);
  }

  return $fc->{_cdna_start};
}

sub _exon_cdna_end {
  my ($self, $exon, $tr) = @_;

  my $tr_stable_id = $tr->stable_id;
  my $fc = $exon->{_variation_effect_feature_cache}->{$tr_stable_id} ||= {};

  if(!exists($fc->{_cdna_end})) {
    $fc->{_cdna_end} = $exon->cdna_end($tr);
  }

  return $fc->{_cdna_end};
}

# same for $transcript->feature_Slice
# need to be careful here in case the transcript has moved slice
# you never know!
sub _transcript_feature_Slice {
  my ($self, $tr) = @_;

  my $fc = $tr->{_variation_effect_feature_cache} ||= {};

  # check that we haven't moved slice
  my $curr_slice_ref = sprintf('%s', $tr->slice());
  my $prev_slice_ref = $fc->{slice_ref};

  if(
    !exists($fc->{feature_Slice}) ||
    $fc->{slice_ref} && $fc->{slice_ref} ne $curr_slice_ref
  ) {

    # log the reference of this slice
    $fc->{slice_ref} = $curr_slice_ref;
    $fc->{feature_Slice} = $tr->feature_Slice();
  }

  return $fc->{feature_Slice};
}


1;


