=head1 LICENSE

 Copyright (c) 1999-2012 The European Bioinformatics Institute and
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
    print "amino acid change: ", $tva->peptide_allele_string, "\n";
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
use Bio::EnsEMBL::Utils::Exception qw(throw warning deprecate);
use Bio::EnsEMBL::Variation::Utils::Sequence qw(hgvs_variant_notation format_hgvs_string);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(within_cds within_intron stop_lost affects_start_codon);

use base qw(Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele Bio::EnsEMBL::Variation::BaseTranscriptVariationAllele);


our $DEBUG = 0;

sub new_fast {
    my ($class, $hashref) = @_;
    
    # swap a transcript_variation argument for a variation_feature_overlap one

    if ($hashref->{transcript_variation}) {
        $hashref->{variation_feature_overlap} = delete $hashref->{transcript_variation};
    }
    
    # and call the superclass

    return $class->SUPER::new_fast($hashref);
}

=head2 transcript_variation

  Description: Get/set the associated TranscriptVariation
  Returntype : Bio::EnsEMBL::Variation::TranscriptVariation
  Exceptions : throws if the argument is the wrong type
  Status     : At Risk

=cut

sub transcript_variation {
    my ($self, $tv) = @_;
    assert_ref($tv, 'Bio::EnsEMBL::Variation::TranscriptVariation') if $tv;
    return $self->variation_feature_overlap($tv);
}

=head2 variation_feature

  Description: Get the associated VariationFeature
  Returntype : Bio::EnsEMBL::Variation::VariationFeature
  Exceptions : none
  Status     : At Risk

=cut

sub variation_feature {
    my $self = shift;
    return $self->transcript_variation->variation_feature;
}

=head2 pep_allele_string

  Description: Return a '/' delimited string of the reference peptide and the 
               peptide resulting from this allele, or a single peptide if this
               allele does not change the peptide (e.g. because it is synonymous)
  Returntype : string or undef if this allele is not in the CDS
  Exceptions : none
  Status     : At Risk

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
  Status     : At Risk

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
    
    my $ref_display_codon = $self->transcript_variation->get_reference_TranscriptVariationAllele->display_codon;
    
    return undef unless $ref_display_codon;
    
    return $ref_display_codon.'/'.$display_codon;
}

=head2 peptide

  Description: Return the amino acid sequence that this allele is predicted to result in
  Returntype : string or undef if this allele is not in the CDS or is a frameshift
  Exceptions : none
  Status     : At Risk

=cut

sub peptide {
    my ($self, $peptide) = @_;
    
    $self->{peptide} = $peptide if $peptide;
    
    unless ($self->{peptide}) {

        return undef unless $self->seq_is_unambiguous_dna;
        
        if (my $codon = $self->codon) {
            
            # the codon method can set the peptide in some circumstances 
            # so check here before we try an (expensive) translation
            return $self->{peptide} if $self->{peptide};
            
            # translate the codon sequence to establish the peptide allele
            
            # allow for partial codons - split the sequence into whole and partial
            # e.g. AAAGG split into AAA and GG            
            my $whole_codon   = substr($codon, 0, int(length($codon) / 3) * 3);
            my $partial_codon = substr($codon, int(length($codon) / 3) * 3);
            
            my $pep = '';
            
            if($whole_codon) {
                # for mithocondrial dna we need to to use a different codon table
                my $codon_table = $self->transcript_variation->_codon_table;
                
                my $codon_seq = Bio::Seq->new(
                    -seq        => $whole_codon,
                    -moltype    => 'dna',
                    -alphabet   => 'dna',
                );
                
                $pep .= $codon_seq->translate(undef, undef, undef, $codon_table)->seq;
            }
            
            if($partial_codon) {
                $pep .= 'X';
            }
            
            $pep ||= '-';
           
            $self->{peptide} = $pep;
        }
    }
    
    return $self->{peptide};
}

=head2 codon

  Description: Return the codon sequence that this allele is predicted to result in
  Returntype : string or undef if this allele is not in the CDS or is a frameshift
  Exceptions : none
  Status     : At Risk

=cut

sub codon {
    my ($self, $codon) = @_;
    
    $self->{codon} = $codon if defined $codon;
    
    my $tv = $self->transcript_variation;      
    
    return undef unless $tv->translation_start && $tv->translation_end;
   
    return undef unless $self->seq_is_dna;
    
    unless ($self->{codon}) {
      
        # try to calculate the codon sequence
    
        my $seq = $self->feature_seq;
        
        $seq = '' if $seq eq '-';
        
        # calculate necessary coords and lengths
        
        my $codon_cds_start = $tv->translation_start * 3 - 2;
        my $codon_cds_end   = $tv->translation_end * 3;
        my $codon_len       = $codon_cds_end - $codon_cds_start + 1;
        my $vf_nt_len       = $tv->cds_end - $tv->cds_start + 1;
        my $allele_len      = $self->seq_length;
        
        if ($allele_len != $vf_nt_len) {
            if (abs($allele_len - $vf_nt_len) % 3) {
                # this is a frameshift variation, we don't attempt to 
                # calculate the resulting codon or peptide change as this 
                # could get quite complicated 
                return undef;
            }
        }

        # splice the allele sequence into the CDS
        
        my $cds = $tv->_translateable_seq;
    
        substr($cds, $tv->cds_start-1, $vf_nt_len) = $seq;
        
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
  Status     : At Risk

=cut

sub display_codon {
    my $self = shift;

    unless ($self->{_display_codon}) {

        if ($self->codon && defined $self->transcript_variation->codon_position) {
            
            my $display_codon = lc $self->codon;

            # if this allele is an indel then just return all lowercase
            
            if ($self->feature_seq ne '-') {
                
                # codon_position is 1-based, while substr assumes the string starts at 0
                
                my $pos = $self->transcript_variation->codon_position - 1;

                my $len = length $self->feature_seq;

                substr($display_codon, $pos, $len) = uc substr($display_codon, $pos, $len);
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
    
    unless ($self->{$analysis}->{prediction}) {
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

    unless ($self->{$analysis}->{score}) {
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
    
    unless ($self->{sift_prediction}) {
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

    unless ($self->{sift_score}) {
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
            
            my ($prediction, $score) = $matrix->get_prediction(
                $self->transcript_variation->translation_start,
                $self->peptide,
            );

            return wantarray ? ($prediction, $score) : $prediction;
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

=head2 hgvs_coding

  Description: Return a string representing the CDS-level effect of this allele in HGVS format
  Returntype : string or undef if this allele is not in the CDS 
  Exceptions : none
  Status     : At Risk

=cut

sub hgvs_coding {

    deprecate('HGVS coding support has been moved to hgvs_transcript. This method will be removed in the next release.');
    return hgvs_transcript(@_);
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
    return $self->{hgvs_transcript}        if defined $self->{hgvs_transcript} ;

    my $variation_feature_sequence  = $self->variation_feature_seq() ;

    ### don't try to handle odd characters
    return undef if $variation_feature_sequence =~ m/[^ACGT\-]/ig;

    ### no result for reference allele
    return undef if $self->is_reference ==1; 

    ### else evaluate
       
    ### get reference sequence strand
    my $refseq_strand = $self->transcript_variation->transcript->strand();

    my $var_name = $self->transcript_variation->variation_feature->variation_name();

    if($DEBUG ==1){    
    print "\nHGVS transcript: Checking $var_name refseq strand => $refseq_strand seq name : " . $self->transcript_variation->transcript_stable_id() . " var strand " . $self->transcript_variation->variation_feature->strand() . " vf st " . $self->variation_feature->strand()  ." seqname: " . $self->variation_feature->seqname()  ." seq: " . $self->variation_feature_seq ."\n";
    }

    my $hgvs_notation ; ### store components of HGVS string in hash
   
    ### vf strand is relative to transcript or transcript slice
    if( $self->transcript_variation->variation_feature->strand() <0 && $refseq_strand >0 ||
    $self->transcript_variation->variation_feature->strand() >0 && $refseq_strand < 0
    ){    
    reverse_comp(\$variation_feature_sequence);
    }
        
    ### need to get ref seq from slice transcript is on for intron labelling    
    my ($slice_start, $slice_end, $slice) = $self->_var2transcript_slice_coords();
    
   
    unless($slice_start){
    #warn "TVA: VF not within transcript - no HGVS\n";
    return undef;
    }
        
    ### decide event type from HGVS nomenclature   
    $hgvs_notation = hgvs_variant_notation(  $variation_feature_sequence,    ### alt_allele,
                         $slice->seq(),                  ### using this to extract ref allele
                         $slice_start,
                         $slice_end,
                         "",
                         "",
                         $var_name 
    );
                       
    ### This should not happen
    unless($hgvs_notation->{'type'}){
    #warn "Error - not continuing; no HGVS annotation\n";
    return undef;
    } 
 
    ## compare calculated reference base to input reference base to flag incorrect input
    my $ref_al_for_checking  = $self->transcript_variation->get_reference_TranscriptVariationAllele->variation_feature_seq();
    if( $self->transcript_variation->variation_feature->strand() <0 && $refseq_strand >0 ||
        $self->transcript_variation->variation_feature->strand() >0 && $refseq_strand < 0
       ){    
        reverse_comp(\$ref_al_for_checking);
    }
    $ref_al_for_checking =~ s/-//;
    unless( $hgvs_notation->{ref}  eq $ref_al_for_checking ||
            $hgvs_notation->{type} eq 'dup' ||
            $hgvs_notation->{type} eq 'inv'){
        warn "\nError - calculated reference ($hgvs_notation->{ref}) and input reference ($ref_al_for_checking) disagree - skipping HGVS transcript for $var_name\n";
    }

    ### create reference name - transcript name & seq version
    my $stable_id = $self->transcript_variation->transcript_stable_id();    
    $stable_id .= "." . $self->transcript_variation->transcript->version() unless $stable_id =~ /\.\d+$/;
    $hgvs_notation->{'ref_name'} =  $stable_id;
  

    ### get position relative to transcript features [use HGVS coords not variation feature coords due to dups]
    $hgvs_notation->{start} = $self->_get_cDNA_position( $hgvs_notation->{start} );
    $hgvs_notation->{end}   = $self->_get_cDNA_position( $hgvs_notation->{end} );
    

 # Make sure that start is always less than end
    my ($exon_start_coord, $intron_start_offset) = $hgvs_notation->{start} =~ m/(\-?[0-9]+)\+?(\-?[0-9]+)?/;
    my ($exon_end_coord,   $intron_end_offset)   = $hgvs_notation->{end} =~ m/(\-?[0-9]+)\+?(\-?[0-9]+)?/;
    $intron_start_offset ||= 0;
    $intron_end_offset   ||= 0;

    ($hgvs_notation->{start},$hgvs_notation->{end}) = ($hgvs_notation->{end},$hgvs_notation->{start}) if 
    (($exon_start_coord + $intron_start_offset) > ($exon_end_coord + $intron_end_offset));



    if($self->transcript->cdna_coding_start()){
    $hgvs_notation->{'numbering'} = "c";  ### set 'c' if transcript is coding 
    }    
    else{
    $hgvs_notation->{'numbering'} = "n";  ### set 'n' if transcript non-coding 
    } 
    ### generic formatting 
    $self->{hgvs_transcript} =  format_hgvs_string( $hgvs_notation);

    if($DEBUG ==1){print "HGVS notation for var " . $self->transcript_variation->variation_feature->variation_name() . " from hgvs transcript : " . $self->{hgvs_transcript} . " \n";}
    
    return $self->{hgvs_transcript};
   
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
    
    if($DEBUG ==1){print "\nStarting hgvs_protein with ".  $self->transcript_variation->variation_feature->variation_name() . "\n"; }

    ### set if string supplied
    $self->{hgvs_protein} = $notation  if defined $notation;
    
    ### return if set
    return $self->{hgvs_protein}       if defined $self->{hgvs_protein} ;
    
    ### don't try to handle odd characters
    return undef if $self->variation_feature_seq() =~ m/[^ACGT\-]/ig;

    ### no HGVS annotation for reference allele
    return undef if $self->is_reference();

    ### no HGVS protein annotation for variants outside translated region 
    return undef unless ($self->transcript_variation->translation_start() && $self->transcript_variation->translation_end()); 
         
    
    
    ### get reference sequence [add seq version to transcript name]
    my $stable_id = $self->transcript_variation->transcript->translation->display_id();
    $stable_id .= "." . $self->transcript_variation->transcript->translation->version() unless $stable_id =~ /\.\d+$/;
    $hgvs_notation->{ref_name} =  $stable_id;

    $hgvs_notation->{'numbering'} = 'p';

    ### get default reference location [changed later in some cases eg. duplication]
    $hgvs_notation->{start}   = $self->transcript_variation->translation_start();
    $hgvs_notation->{end}     = $self->transcript_variation->translation_end();  


    ## get default reference & alt peptides  [changed later to hgvs format]
    $hgvs_notation->{alt} = $self->peptide;
    $hgvs_notation->{ref} = $self->transcript_variation->get_reference_TranscriptVariationAllele->peptide;    

    if(defined $hgvs_notation->{alt} && defined $hgvs_notation->{ref}){
        $hgvs_notation = _clip_alleles( $hgvs_notation);
    }

    #### define type - types are different for protein numbering
    $hgvs_notation  = $self->_get_hgvs_protein_type($hgvs_notation);

    ##### Convert ref & alt peptides taking into account HGVS rules
    $hgvs_notation = $self->_get_hgvs_peptides($hgvs_notation);

    ##### String formatting
    $self->{hgvs_protein}  =  $self->_get_hgvs_protein_format($hgvs_notation);
    return $self->{hgvs_protein} ;   
}

### HGVS: format protein string
sub _get_hgvs_protein_format{

    my $self          = shift;
    my $hgvs_notation = shift;

    if ((defined  $hgvs_notation->{ref} && defined $hgvs_notation->{alt} && 
     $hgvs_notation->{ref} eq $hgvs_notation->{alt}) || 
      (defined  $hgvs_notation->{type} && $hgvs_notation->{type} eq "=")){

        ### no protein change - return transcript nomenclature with flag for neutral protein consequence
	if(defined $self->hgvs_transcript()){
	    $hgvs_notation->{'hgvs'} = $self->hgvs_transcript() . "(p.=)";
	    return $hgvs_notation->{'hgvs'} ;
	}
	else{
	    return undef;
	}
    }

    ### all start with refseq name & numbering type
    $hgvs_notation->{'hgvs'} = $hgvs_notation->{'ref_name'} . ":" . $hgvs_notation->{'numbering'} . ".";    

    ### handle stop_lost seperately regardless of cause by del/delins => p.XposAA1extnum_AA_to_stop
    if(stop_lost($self)){  ### if deletion of stop add extX and number of new aa to alt
        $hgvs_notation->{alt} = substr($hgvs_notation->{alt}, 0, 3);
        if($hgvs_notation->{type} eq "del"){
            my $aa_til_stop =  _stop_loss_extra_AA($self,$hgvs_notation->{start}-1, "del");
            if(defined  $aa_til_stop){
                $hgvs_notation->{alt} .=  "extX" . $aa_til_stop;
            }
        }
        elsif($hgvs_notation->{type} eq ">"){
            my $aa_til_stop =  _stop_loss_extra_AA($self,$hgvs_notation->{start}-1, "subs");
            if(defined  $aa_til_stop){  
              $hgvs_notation->{alt} .=  "extX" . $aa_til_stop;
            }
        }
        else{
           # warn "TVA: stop loss for type $hgvs_notation->{type}  not caught \n";
        }
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
    
        #### list first and last aa in reference only
        my $ref_pep_first = substr($hgvs_notation->{ref}, 0, 3);
        my $ref_pep_last;
        if(substr($hgvs_notation->{ref}, -1, 1) eq "X"){
            $ref_pep_last ="X";
        }
        else{
            $ref_pep_last  = substr($hgvs_notation->{ref}, -3, 3);
         }
         if($hgvs_notation->{ref} =~ /X$/){
            ### For stops & add extX & distance to next stop to alt pep
            my $aa_til_stop =  _stop_loss_extra_AA($self,$hgvs_notation->{start}-1, "loss");
            if(defined  $aa_til_stop){  
               $hgvs_notation->{alt} .="extX" . $aa_til_stop;
            }
         }
         if($hgvs_notation->{start} == $hgvs_notation->{end} && $hgvs_notation->{type} eq "delins"){       
             $hgvs_notation->{'hgvs'} .= $ref_pep_first . $hgvs_notation->{start} . $hgvs_notation->{end} . $hgvs_notation->{type} . $hgvs_notation->{alt} ;
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
       if(defined $hgvs_notation->{alt} && $hgvs_notation->{alt} eq "X"){ ## stop gained
           $hgvs_notation->{'hgvs'} .= $hgvs_notation->{ref} . $hgvs_notation->{start}  .  $hgvs_notation->{alt} ;

       }
       else{ ## not immediate stop - count aa until next

          my $aa_til_stop =  _stop_loss_extra_AA($self, $hgvs_notation->{start}-1, "fs");
          if($aa_til_stop ){
              ### use long form if new stop found 
             $hgvs_notation->{'hgvs'} .= $hgvs_notation->{ref} . $hgvs_notation->{start}  .  $hgvs_notation->{alt} . "fsX$aa_til_stop"  ;
          }
          else{
            ### otherwise use short form 
            $hgvs_notation->{'hgvs'} .= $hgvs_notation->{ref} . $hgvs_notation->{start}  . "fs";
          }
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
sub _get_hgvs_protein_type{

    my $self = shift;
    my $hgvs_notation = shift;
   
    ### get allele length
    my ($ref_length, $alt_length ) = $self->_get_allele_length();    


    if( defined $hgvs_notation->{ref} && defined $hgvs_notation->{alt} ){
    ### Run type checks on peptides if available

      if($hgvs_notation->{alt} eq $hgvs_notation->{ref} ){
          #### synonymous indicated by =
          $hgvs_notation->{type} = "=";
      }
      elsif( length($hgvs_notation->{ref}) ==1 && length($hgvs_notation->{alt}) ==1 ) {
 
          $hgvs_notation->{type} = ">";
      }
      elsif($hgvs_notation->{ref} eq "-" || $hgvs_notation->{ref} eq "") {

          $hgvs_notation->{type} = "ins";
      }
      elsif($hgvs_notation->{alt} eq "" ) {

          $hgvs_notation->{type} = "del";
      }
      elsif( (length($hgvs_notation->{alt}) >  length($hgvs_notation->{ref}) ) &&       
           $hgvs_notation->{alt} =~ /$hgvs_notation->{ref}/ ) {
          ### capture duplication event described as TTT/TTTTT 
          $hgvs_notation->{type} = "dup";        
      }

      elsif( (length($hgvs_notation->{alt}) >0 && length($hgvs_notation->{ref}) >0) &&
             (length($hgvs_notation->{alt}) ne length($hgvs_notation->{ref}) ) ) {
          $hgvs_notation->{type} = "delins";
      }
      else{
          $hgvs_notation->{type} = ">";
      }
    }
    
    ## not a simple event - check DNA strings
    elsif($ref_length ne $alt_length && ($ref_length  - $alt_length)%3 !=0 ){
        $hgvs_notation->{type} = "fs";
    }    
    
    elsif($alt_length >1  ){
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
    else{
        #print STDERR "DEBUG ".$self->variation_feature->start."\n";
        #warn "Cannot define protein variant type [$ref_length  - $alt_length]\n";
    }
    return $hgvs_notation ;

}

### HGVS: get reference & alternative peptide 
sub _get_hgvs_peptides{

  my $self          = shift;
  my $hgvs_notation = shift;

  if($hgvs_notation->{type} eq "fs"){
    ### ensembl alt/ref peptides not the same as HGVS alt/ref - look up seperately
    $hgvs_notation = $self->_get_fs_peptides($hgvs_notation);    
  }
  elsif($hgvs_notation->{type} eq "ins" ){

    ### HGVS ref are peptides flanking insertion
    $hgvs_notation->{ref} = $self->_get_surrounding_peptides($hgvs_notation->{start});

    if( $hgvs_notation->{alt} =~/\*/){
        ## inserted bases after stop irrelevant; consider as substitution gaining stop MAINTAIN PREVIOUS BEHAVIOUR FOR NOW
        #$hgvs_notation->{alt} =~ s/\*\w+$/\*/;
    }
    else{
        ### Additional check that inserted bases do not duplicate 3' reference sequence [set to type = dup if so]
        $hgvs_notation = $self->_check_for_peptide_duplication($hgvs_notation);
    }
  }
  

  ### Convert peptide to 3 letter code as used in HGVS
  $hgvs_notation->{ref}  = Bio::SeqUtils->seq3(Bio::PrimarySeq->new(-seq => $hgvs_notation->{ref}, -id => 'ref',  -alphabet => 'protein')) || "";
  if( $hgvs_notation->{alt} eq "-"){
    $hgvs_notation->{alt} = "del";
  }
  else{
    $hgvs_notation->{alt}  = Bio::SeqUtils->seq3(Bio::PrimarySeq->new(-seq => $hgvs_notation->{alt}, -id => 'ref',  -alphabet => 'protein')) || "";
  }

  ### handle special cases
  if(    affects_start_codon($self) ){
    #### handle initiator loss -> probably no translation => alt allele is '?'
    $hgvs_notation->{alt}  = "?";    
    $hgvs_notation->{type} = "";
  }

  elsif(  $hgvs_notation->{type} eq "del"){ 
    #### partial last codon    
    $hgvs_notation->{alt} = "del";
  } 
  elsif($hgvs_notation->{type} eq "fs"){
    ### only quote first ref peptide for frameshift
    $hgvs_notation->{ref} = substr($hgvs_notation->{ref},0,3);
  }

  ### set stop to HGVS prefered "X"
  if(defined $hgvs_notation->{ref}){ $hgvs_notation->{ref} =~ s/Ter|Xaa/X/g; }
  if(defined $hgvs_notation->{alt}){ $hgvs_notation->{alt} =~ s/Ter|Xaa/X/g; }
         
  return ($hgvs_notation);                  
                
}

### HGVS:  remove common peptides from alt and ref strings & shift start/end accordingly 
sub _clip_alleles{
    
  my $hgvs_notation = shift;
    
  my $check_alt   = $hgvs_notation->{alt} ;
  my $check_ref   = $hgvs_notation->{ref} ;
  my $check_start = $hgvs_notation->{start};
    
  ### strip same bases from start of string
  for (my $p =0; $p <length ($hgvs_notation->{ref}); $p++){
  my $check_next_ref = substr( $check_ref, 0, 1);
  my $check_next_alt = substr( $check_alt, 0, 1);
    
    if($check_next_ref eq  "*" && $check_next_alt eq "*"){
        ### stop re-created by variant - no protein change
        $hgvs_notation->{type} = "=";

        return($hgvs_notation);
    }
    
    if($check_next_ref eq  $check_next_alt){
        $check_start++;
        $check_ref  = substr( $check_ref, 1);
        $check_alt  = substr( $check_alt, 1);
    }
    else{
        last;
    }
  }
  #### strip same bases from end of string
  for (my $q =0; $q < length ($check_ref); $q++){
    my $check_next_ref = substr( $check_ref, -1, 1);
    my $check_next_alt = substr( $check_alt, -1, 1);
    if($check_next_ref eq  $check_next_alt){
        chop $check_ref;
        chop $check_alt;
    }
    else{
        last;
    }
  }

  $hgvs_notation->{alt}   = $check_alt;
  $hgvs_notation->{ref}   = $check_ref;

  ### check if clipping force type change & adjust location
  if(length ($check_ref) == 0  && length ($check_alt) >= 1){
    ### re-set as ins not delins  
    $hgvs_notation->{type} ="ins";
        ### insertion between ref base and next => adjust next         
    if($hgvs_notation->{end} == $hgvs_notation->{start} ){$hgvs_notation->{end} = $check_start;    }
#    $hgvs_notation->{start} = $check_start;    
  }
  elsif(length ($check_ref) >=1  && length ($check_alt) == 0){
    ### re-set as del not delins  
    $hgvs_notation->{type}  ="del";        
    $hgvs_notation->{start} = $check_start;
  }
  else{
    #### save trimmed peptide string & increased start position
    $hgvs_notation->{start} = $check_start;

  }        

    return $hgvs_notation;
}
    




#### HGVS: check allele lengths to look for frameshifts
sub _get_allele_length{

  my $self       = shift;
  my $ref_length = 0;
  my $alt_length = 0;

  my $al_string = $self->allele_string();
  my $ref_allele   = (split/\//, $al_string)[0];
  $ref_allele =~ s/\-//;
  $ref_length    =  length $ref_allele;

  my $alt_allele = $self->variation_feature_seq();
  $alt_allele =~ s/\-//;
  $alt_length    =  length  $alt_allele;

  return ($ref_length, $alt_length );
  
}

### HGVS: list first different peptide [may not be first changed codon]
sub _get_fs_peptides{

  my $self    = shift;
  my $hgvs_notation = shift;

  ### get CDS with alt variant
  my $alt_cds = $self->_get_alternate_cds();
  return undef unless defined($alt_cds);

  #### get new translation
  my $alt_trans = $alt_cds->translate()->seq();
  
  ### get changed end (currently in single letter AA coding)    
  my $ref_trans  = $self->transcript->translate()->seq();
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
          ### variation at stop codon, but maintains stop codon
          return ($hgvs_notation);
      }
      last if $hgvs_notation->{ref} ne $hgvs_notation->{alt};
      $hgvs_notation->{start}++;
    }
    return ($hgvs_notation);

}

#### HGVS: if variant is an insertion, ref pep is initially "-", so seek AA before and after insertion
sub _get_surrounding_peptides{

  my $self    = shift;
  my $ref_pos = shift; 
  my $ref_trans  = $self->transcript->translate()->seq();

  my $end = substr($ref_trans, $ref_pos-1);
  my $ref_string = substr($ref_trans, $ref_pos-1, 2);

  return ($ref_string);

}


#### HGVS: alternate CDS needed to check for new stop when variant disrupts 'reference' stop
sub _get_alternate_cds{
    
  my $self = shift;

  ### get reference sequence
  my $reference_cds_seq = $self->transcript->translateable_seq();
  
  return undef unless defined($self->transcript_variation->cds_start) && defined($self->transcript_variation->cds_end());

  ### get sequences upstream and downstream of variant
  my $upstream_seq   =  substr($reference_cds_seq, 0, ($self->transcript_variation->cds_start() -1) );
  my $downstream_seq =  substr($reference_cds_seq, ($self->transcript_variation->cds_end() ) );

  ### fix alternate allele if deletion or on opposite strand
  my $alt_allele  = $self->variation_feature_seq();
  $alt_allele  =~ s/\-//;
  if( $self->transcript_variation->variation_feature->strand() <0 && $self->transcript_variation->transcript->strand() >0 ||
      $self->transcript_variation->variation_feature->strand() >0 && $self->transcript_variation->transcript->strand() < 0
      ){    
      reverse_comp(\$alt_allele) ;
  }
  ### build alternate seq
  my $alternate_seq  = $upstream_seq . $alt_allele . $downstream_seq ;

  ### create seq obj with alternative allele in the CDS sequence
  my $alt_cds =Bio::PrimarySeq->new(-seq => $alternate_seq,  -id => 'alt_cds', -alphabet => 'dna');

  ### append UTR if available as stop may be disrupted
  my $utr = $self->transcript_variation->transcript->three_prime_utr();

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
sub _check_for_peptide_duplication{
    
   my $self = shift;
   my $hgvs_notation= shift;

   ##### get reference sequence
   my $reference_cds_seq = $self->transcript->translateable_seq();

   ##### get sequence upstream of variant
   my $upstream_seq   =  substr($reference_cds_seq, 0, ($self->transcript_variation->cds_start() -1) );
   
   ##### create translation to check against inserted peptides 
   my $upstream_cds =Bio::PrimarySeq->new(-seq => $upstream_seq,  -id => 'alt_cds', -alphabet => 'dna');
   my $upstream_trans = $upstream_cds->translate()->seq();
   
   ## Test whether alt peptide matches the reference sequence just before the variant
   my $test_new_start = $hgvs_notation->{'start'} - length($hgvs_notation->{'alt'}) -1 ;
   my $test_seq       =  substr($upstream_trans, $test_new_start, length($hgvs_notation->{'alt'}));

   if ( $test_new_start >= 0 && $test_seq eq $hgvs_notation->{alt}) {

	  $hgvs_notation->{type}   = 'dup';
	  $hgvs_notation->{end}    = $hgvs_notation->{start} -1;
	  $hgvs_notation->{start} -= length($hgvs_notation->{alt});

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
  
  my $ref_temp  =  $self->transcript->translate()->seq;
  my $ref_len = length($ref_temp);
  
  if($DEBUG==1){ 
    print "alt translated:\n" . $alt_trans->seq() . "\n";
    print "ref translated:\n$ref_temp\n";;
  }
  
  #### Find the number of residues that are translated until a termination codon is encountered
  if ($alt_trans->seq() =~ m/\*/) {
    if($DEBUG==1){print "Got $+[0] aa before stop, var event at $ref_var_pos \n";}
  
    if($test eq "fs" ){
      ### frame shift - count from first AA effected by variant to stop
      $extra_aa = $+[0] - $ref_var_pos;
      if($DEBUG==1){ print "Stop change ($test): found $extra_aa amino acids before fs stop [ $+[0] - peptide ref_start: $ref_var_pos )]\n";}
    }
  
    else{
      $extra_aa = $+[0]  - 1 - $ref_len;
      if($DEBUG==1){ print "Stop change ($test): found $extra_aa amino acids before next stop [ $+[0] - 1 -normal stop $ref_len)]\n";}        
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
    
    #ÊThe rna and mitochondrial modes have not yet been implemented, so return undef in case we get a call to these
    return undef if ($reference =~ m/rna|mitochondrial/);
    
    my $sub = qq{hgvs_$reference};
    
    $self->{$sub} = $notation if defined $notation;
    
    unless ($self->{$sub}) {

    # Use the transcript this VF is on as the reference feature
        my $reference_feature = $self->transcript;
        # If we want genomic coordinates, the reference_feature should actually be the slice for the underlying seq_region
        $reference_feature = $reference_feature->slice->seq_region_Slice if ($reference eq 'genomic');

        # Calculate the HGVS notation on-the-fly and pass it to the TranscriptVariation in order to distribute the result to the other alleles
        $self->transcript_variation->$sub($self->variation_feature->get_all_hgvs_notations($reference_feature,substr($reference,0,1),undef,undef,$self->transcript_variation));
    }
    
    return $self->{$sub};
}


### HGVS: move variant to transcript slice
sub _var2transcript_slice_coords{

  my $self = shift;
  
  my $ref_slice = $self->transcript->feature_Slice();  #### returns with strand same as feature
  my $tr_vf     = $self->variation_feature->transfer($ref_slice);
          
  # Return undef if this VariationFeature does not fall within the supplied feature.
  return undef if (!defined $tr_vf ||
                   $tr_vf->start  < 1 || 
                   $tr_vf->end    < 1 || 
                   $tr_vf->start  > ($self->transcript->end - $self->transcript->start + 1) || 
                   $tr_vf->end    > ($self->transcript->end - $self->transcript->start + 1)); 
  
  return( $tr_vf->start() , $tr_vf->end(), $ref_slice);
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
    
    my $transcript = $self->transcript();
    my $strand     = $transcript->strand();    

    #### TranscriptVariation start/stop coord relative to transcript 
    #### Switch to chromosome coordinates taking into account strand
    $position = ( $strand > 0 ? 
          ( $self->transcript->start() + $position - 1 )  :   
          ( $self->transcript->end()   - $position + 1));


            
    # Get all exons and sort them in positional order
    my @exons   = sort {$a->start() <=> $b->start()} @{$transcript->get_all_Exons()};
    my $n_exons = scalar(@exons);

    my $cdna_position;
  # Loop over the exons and get the coordinates of the variation in exon+intron notation
    for (my $i=0; $i<$n_exons; $i++) {
    
      # Skip if the start point is beyond this exon
      next if ($position > $exons[$i]->end());
    

      # EXONIC:  If the start coordinate is within this exon
      if ($position >= $exons[$i]->start()) {
      # Get the cDNA start coordinate of the exon and add the number of nucleotides from the exon boundary to the variation
      # If the transcript is in the opposite direction, count from the end instead
      $cdna_position = $exons[$i]->cdna_start($transcript) + ( $strand > 0 ? 
                                   ( $position - $exons[$i]->start ) : 
                                   ( $exons[$i]->end() - $position ) 
                                   );
      last;  #### last exon checked
    }
      ## INTRONIC
      # Else the start coordinate is between this exon and the previous one, determine which one is closest and get coordinates relative to that one
      else {

      my $updist   = ($position - $exons[$i-1]->end());
      my $downdist = ($exons[$i]->start() - $position);

      # If the distance to the upstream exon is the shortest, or equal and in the positive orientation, use that
      if ($updist < $downdist || ($updist == $downdist && $strand >= 0)) {
          
          # If the orientation is reversed, we should use the cDNA start and a '-' offset
          $cdna_position = ($strand >= 0 ? 
                $exons[$i-1]->cdna_end($transcript)   . '+' : 
                $exons[$i-1]->cdna_start($transcript) . '-') 
          . $updist;
      }
      # Else if downstream is shortest...
      else {
          # If the orientation is reversed, we should use the cDNA end and a '+' offset
          $cdna_position = ($strand >= 0 ? 
                $exons[$i]->cdna_start($transcript) . '-' : 
                $exons[$i]->cdna_end($transcript)   . '+') 
          . $downdist;
      }
      last; ## last exon checked
      }
  }
  
  # Shift the position to make it relative to the start & stop codons
  my $start_codon  =  $transcript->cdna_coding_start();
  my $stop_codon   =  $transcript->cdna_coding_end();
  
  # Disassemble the cDNA coordinate into the exon and intron parts
    ### just built this now taking it appart again
  my ($cdna_coord, $intron_offset) = $cdna_position =~ m/([0-9]+)([\+\-][0-9]+)?/;
  

  # Start by correcting for the stop codon
  if (defined($stop_codon) && $cdna_coord > $stop_codon) {
    # Get the offset from the stop codon
    $cdna_coord -= $stop_codon;
    # Prepend a * to indicate the position is in the 3' UTR
    $cdna_coord = '*' . $cdna_coord;
  }
  elsif (defined($start_codon)) {
    # If the position is beyond the start codon, add 1 to get the correct offset
    $cdna_coord += ($cdna_coord >= $start_codon);
    # Subtract the position of the start codon
    $cdna_coord -= $start_codon;
  }

  # Re-assemble the cDNA position  [ return exon num & offset & direction for intron eg. 142+363] 
  $cdna_position = $cdna_coord . (defined($intron_offset) ? $intron_offset : '');
  return $cdna_position;
}


1;


