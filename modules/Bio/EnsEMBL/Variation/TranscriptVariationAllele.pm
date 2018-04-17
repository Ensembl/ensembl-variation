=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

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
  Status     : Stable

=cut

sub display_codon_allele_string {
    my ($self) = @_;
    
    my $display_codon = $self->display_codon;
    
    return undef unless $display_codon;
    
    my $ref_display_codon = $self->transcript_variation->get_reference_TranscriptVariationAllele->display_codon;
    
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

    my ($tv_tr_start, $tv_tr_end) = ($tv->translation_start, $tv->translation_end);

    unless($tv_tr_start && $tv_tr_end && $self->seq_is_dna) {
      return $self->{codon};
    }
  
    # try to calculate the codon sequence
    my $seq = $self->feature_seq;
    
    $seq = '' if $seq eq '-';
    
    # calculate necessary coords and lengths
    
    my $codon_cds_start = $tv_tr_start * 3 - 2;
    my $codon_cds_end   = $tv_tr_end * 3;
    my $codon_len       = $codon_cds_end - $codon_cds_start + 1;
    my $vf_nt_len       = $tv->cds_end - $tv->cds_start + 1;
    my $allele_len      = $self->seq_length;
    
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
      $cds = $cds_obj->seq();
    }

    else {
      # splice the allele sequence into the CDS
      $cds = $tv->_translateable_seq;
    
      substr($cds, $tv->cds_start-1, $vf_nt_len) = $seq;
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
               if this is a non-synonymous mutation and a prediction is available, undef
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
  Returntype : float between 0 and 1 if this is a non-synonymous mutation and a prediction is 
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
  Returntype : string (one of 'tolerated', 'deleterious') if this is a non-synonymous 
               mutation and a prediction is available, undef otherwise
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
  Returntype : float between 0 and 1 if this is a non-synonymous mutation and a prediction is 
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

  ## create new transcript variation object as position may be different
  ## for most variants (SNPs etc) this will actually just return $self
  ## logic is neater this way though
  my $hgvs_tva = $self->_hgvs_tva($tr, $tv, $vf);

  ## return if a new transcript_variation_allele is not available - variation outside transcript
  return undef unless defined $hgvs_tva && defined $hgvs_tva->base_variation_feature_overlap;

  unless (defined  $self->{_slice_start} ){
  	print "Exiting hgvs_transcript: no slice start position for $var_name in trans" . $tr_stable_id . "\n " if $DEBUG == 1 ;
  	return undef;
  }

  ## this may be different to the input one for insertions/deletions
  my $variation_feature_sequence = $hgvs_tva->variation_feature_seq();
  if($variation_feature_sequence && $vf->strand() != $refseq_strand) {    
    reverse_comp(\$variation_feature_sequence) ;
  };

  ### decide event type from HGVS nomenclature   
  print "sending alt: $variation_feature_sequence &  $self->{_slice_start} -> $self->{_slice_end} for formatting\n" if $DEBUG ==1;
  $hgvs_notation = hgvs_variant_notation(
    $variation_feature_sequence,    ### alt_allele,
    $self->{_slice}->seq(),                             ### using this to extract ref allele
    $self->{_slice_start},
    $self->{_slice_end},
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
  $hgvs_notation->{start} = $hgvs_tva->_get_cDNA_position( $hgvs_notation->{start} );
  $hgvs_notation->{end}   = $same_pos ? $hgvs_notation->{start} : $hgvs_tva->_get_cDNA_position( $hgvs_notation->{end} );
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

=head2 hgvs_transcript_reference

  Description: Return a string representing the reference sequence as could be used in HGVS notation
               Useful for deletions where the recommended format formerly used the reference sequence
               but no longer does.
  Returntype : string or undef if no reference allele is available.
  Exceptions : none
  Status     : At Risk

=cut

sub hgvs_transcript_reference{

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
  ## create new tva so position can be changed without impacting ensembl consequence
  my $hgvs_tva = $self->_hgvs_tva();

  ## return if a new transcript_variation_allele is not available - variation outside transcript
  return undef unless defined $hgvs_tva;
  
  my $hgvs_tva_tv = $hgvs_tva->base_variation_feature_overlap;
  return undef unless defined $hgvs_tva_tv;

  my $hgvs_tva_vf = $hgvs_tva_tv->base_variation_feature;
  my $tr          = $hgvs_tva_tv->transcript;
  my $pre         = $hgvs_tva->_pre_consequence_predicates;

  ### no HGVS protein annotation for variants outside translated region 
  unless (
    $pre->{coding} &&
    $hgvs_tva_tv->translation_start() && 
    $hgvs_tva_tv->translation_end()
  ){
    print "Exiting hgvs_protein - variant " . $hgvs_tva_vf->variation_name() . "not within translation\n"  if $DEBUG == 1;
    return undef;
  }
       
  print "proceeding with hgvs prot\n" if $DEBUG == 1;
  print "Checking translation start: " . $hgvs_tva_tv->translation_start() ."\n" if $DEBUG == 1;

  ## checks complete - start building term

  ### get reference sequence and add seq version unless LRG
  $hgvs_notation->{ref_name} = $tr->translation->display_id();
  $hgvs_notation->{ref_name} .= "." . $tr->translation->version() 
    unless ($hgvs_notation->{ref_name}=~ /\.\d+$/ || $hgvs_notation->{ref_name} =~ /LRG/);

  $hgvs_notation->{'numbering'} = 'p';

  ### get default reference location [changed later in some cases eg. duplication]
  $hgvs_notation->{start}   = $hgvs_tva_tv->translation_start();
  $hgvs_notation->{end}     = $hgvs_tva_tv->translation_end();  


  ## get default reference & alt peptides  [changed later to hgvs format]
  $hgvs_notation->{alt} = $hgvs_tva->peptide;
  $hgvs_notation->{ref} = $hgvs_tva_tv->get_reference_TranscriptVariationAllele->peptide;    
  print "Got protein peps: $hgvs_notation->{ref} =>  $hgvs_notation->{alt} (" . $hgvs_tva->codon() .")\n" if $DEBUG ==1;

  
  if(defined $hgvs_notation->{alt} && defined $hgvs_notation->{ref} &&
    ($hgvs_notation->{alt} ne  $hgvs_notation->{ref})){
    $hgvs_notation = _clip_alleles( $hgvs_notation);
  }


  #### define type - types are different for protein numbering
  $hgvs_notation  = $hgvs_tva->_get_hgvs_protein_type($hgvs_notation);
  return undef unless defined $hgvs_notation->{type}; 

  ##### Convert ref & alt peptides taking into account HGVS rules
  $hgvs_notation = $hgvs_tva->_get_hgvs_peptides($hgvs_notation);

  unless($hgvs_notation) {
    $self->{hgvs_protein} = undef;
    return undef;
  }

  ##### String formatting
  return $hgvs_tva->_get_hgvs_protein_format($hgvs_notation);
}

sub _hgvs_tva {
  my ($self, $tr, $tv, $vf) = @_;

  # return self if we've worked out we didn't need to create a new one
  return $self if $self->{_hgvs_tva_is_self};

  if(!exists($self->{_hgvs_tva})) {

    $tv ||= $self->base_variation_feature_overlap;
    $vf ||= $tv->base_variation_feature;
    $tr ||= $tv->transcript;

    my $offset = 0;

    ## need to get ref seq from transcript slice
    my ($slice_start, $slice_end, $slice ) = $self->_var2transcript_slice_coords($tr, $tv, $vf);

    ## no annotation possible if variant outside transcript
    unless ($slice_start) {
      print "couldn't get slice pos for " .$vf->variation_name() . " " . $tr->display_id() . "\n" if $DEBUG ==1; 
      $self->{_hgvs_tva} = undef;
      return $self->{_hgvs_tva};
    }
    
    ## for readability
    my $var_class  =  $vf->var_class();
    $var_class  =~ s/somatic_//;

    ##  only check insertions & deletions & don't move beyond transcript
    if(
      ($var_class eq 'deletion' || $var_class eq 'insertion' ) &&
      $slice_start != $slice->length() &&
      (
        (
          defined $tv->adaptor() && UNIVERSAL::can($tv->adaptor, 'isa') && $tv->adaptor->db ? 
          $tv->adaptor->db->shift_hgvs_variants_3prime()  == 1 :
          $Bio::EnsEMBL::Variation::DBSQL::TranscriptVariationAdaptor::DEFAULT_SHIFT_HGVS_VARIANTS_3PRIME == 1
        )
      )
    ) {

      print "checking position for $var_class, transcript strand is : ".  $tr->strand() ."\n" if $DEBUG ==1;

      my $ref_allele =  $tv->get_reference_TranscriptVariationAllele->variation_feature_seq();
      my $alt_allele =  $self->variation_feature_seq();    

      my $seq_to_check;
      ## sequence to compare is the reference allele for deletion
      $seq_to_check = $ref_allele 	   if $var_class eq 'deletion' ;

      ## sequence to compare is the alt allele
      $seq_to_check = $alt_allele	   if $var_class eq 'insertion';


      my $allele_flipped = 0;
      ## switch allele to transcript strand if necessary
      if(
        ($vf->strand() < 0 && $tr->strand() > 0) ||
        ($vf->strand() > 0 && $tr->strand() < 0)
      ) {
        reverse_comp(\$seq_to_check);
        $allele_flipped = 1;
      }

     	## get reference sequence 3' of variant (slice in transcript strand orientation
      my $from = $slice_end;	    
      my $post_seq = substr($slice->seq(), $from); 
      my $checking_seq = substr( $post_seq , 0,20 );
      print "getting seq: $from to transcript end ($checking_seq )\n" if $DEBUG ==1;

      ## run check
      my $allele;
      ( $allele, $offset ) = get_3prime_seq_offset($seq_to_check, $post_seq );
      print "got allele & offset  $allele, $offset length checked:". length($post_seq) ."\n" if $DEBUG ==1;
      ## correct allele for strand
      reverse_comp(\$allele)   if $allele_flipped == 1;

      if($offset == length($post_seq) ){ ## moved beyond transcript - no annotation	    
        $self->{_hgvs_tva} = undef;
        return $self->{_hgvs_tva};
      }
      elsif($offset > 0 ){
        ## need to move - create new objects with correct position & alleles
        print "moving tva ($allele) by $offset for $var_class, seq checked was length " . length($post_seq) . "\n" if $DEBUG ==1;	  
        print "adding offset $offset to slice pos: $slice_start &  $slice_end\n" if $DEBUG ==1;
        ## increment these & save later - offset positive as slice in transcript orientation
        $slice_start  =  $slice_start + $offset;
        $slice_end    =  $slice_end   + $offset;

        ##  define alleles for variation_feature
        $ref_allele = "-"       if $var_class eq 'insertion';
        $alt_allele = $allele   if $var_class eq 'insertion';	    
        $ref_allele = $allele   if $var_class eq 'deletion';	    
        $alt_allele = "-"       if $var_class eq 'deletion';

        ## make offset negative if rev strand
        $offset = 0 - $offset    if   $tr->strand() <0 ;
        print "sending $ref_allele/$alt_allele & offset:$offset to make_hgvs for $var_class\n" if $DEBUG == 1;

        $self->{_hgvs_tva} = $self->_make_hgvs_tva($ref_allele, $alt_allele, $offset)
      }
      else{
        print "leaving position as input for $var_class after checking \n" if $DEBUG ==1;
        $self->{_hgvs_tva_is_self} = 1;
      }
    }
    
    else {
      print "leaving position as input for $var_class  \n" if $DEBUG ==1;
      $self->{_hgvs_tva_is_self} = 1;
    }

    ## save this to be able to report when HGVS is shifted 
    $self->{_hgvs_offset} = $offset;

    ## add cache of seq/ pos required by c, n and p 
    $self->{_slice_start} = $slice_start;
    $self->{_slice_end}   = $slice_end;
    $self->{_slice}       = $slice;
  }

  return $self->{_hgvs_tva_is_self} ? $self : $self->{_hgvs_tva};
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
  return $self->{_hgvs_offset};
}


sub _make_hgvs_tva {
  my $self       = shift;
  my $ref_allele = shift;
  my $alt_allele = shift;
  my $offset     = shift;

  my $allele_string =  $ref_allele . "/" . $alt_allele;
  
  my $tv        = $self->transcript_variation;
  my $vf        = $tv->variation_feature;

  my $start     = $vf->start() + $offset;
  my $end       = $vf->end() + $offset;

  print "Starting make hgvs tva - vf at $start - $end  $allele_string\n" if $DEBUG ==1;
  print "previous pos :".  $vf->start() ."-" . $vf->end() ."\n" if $DEBUG ==1;

  my $moved_vf =  Bio::EnsEMBL::Variation::VariationFeature->new_fast({
    start          => $start,
    end            => $end,
    allele_string  => $allele_string,
    strand         => $vf->strand(),
    map_weight     => $vf->map_weight(),
    #	-adaptor        => $self->transcript_variation->adaptor->db->get_VariationFeatureAdaptor(),
    variation_name => $vf->variation_name(),
    # variation      => $vf->variation(), ## dont think we need variation, this was running a DB lookup!!!
    slice          => $vf->slice()
  });

  my $transcript = $self->transcript_variation->transcript();
  my $moved_tv = Bio::EnsEMBL::Variation::TranscriptVariation->new(
    -transcript        => $transcript,
    -variation_feature => $moved_vf,
    -no_ref_check      => 1,
    -no_transfer       => 1
  );

  my $hgvs_tva = Bio::EnsEMBL::Variation::TranscriptVariationAllele->new_fast(
    {
     transcript_variation   => $moved_tv,
     variation_feature_seq  => $alt_allele,
     is_reference           => 0
    },
    1									
  );

  return $hgvs_tva;
}

### HGVS: format protein string - runs on $self->{hgvs_tva}
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
    $hgvs_notation->{'hgvs'} .=  $hgvs_notation->{ref} . $hgvs_notation->{start} .  $hgvs_notation->{alt} ;
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

    $hgvs_notation->{alt} = "Ter" if $hgvs_notation->{alt} =~ /^Ter/;

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
sub _get_alternate_cds{
    
  my $self = shift;

  ### get reference sequence
  my $reference_cds_seq = $self->transcript_variation->_translateable_seq();
  
  my $tv = $self->transcript_variation;
  my $vf = $tv->variation_feature;
  my $tr = $tv->transcript;
  
  return undef unless defined($tv->cds_start) && defined($tv->cds_end());

  ### get sequences upstream and downstream of variant
  my $upstream_seq   =  substr($reference_cds_seq, 0, ($tv->cds_start() -1) );
  my $downstream_seq =  substr($reference_cds_seq, ($tv->cds_end() ) );

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
      ($vf_start, $vf_end) = map {($tr_end - $_) + 1} ($vf->end, $vf->start);
    }
    else {
      ($vf_start, $vf_end) = map {($_ - $tr_start) + 1} ($vf->start, $vf->end);
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


