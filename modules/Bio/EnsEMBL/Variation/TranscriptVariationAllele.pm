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

use Bio::EnsEMBL::Variation::Utils::Condel qw(get_condel_prediction);
use Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix qw($AA_LOOKUP);

use base qw(Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele);

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

=head2 transcript

  Description: Get the associated Transcript
  Returntype : Bio::EnsEMBL::Transcript
  Exceptions : none
  Status     : At Risk

=cut

sub transcript {
    my $self = shift;
    return $self->transcript_variation->transcript;
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
    
    return undef unless $tv->translation_start;
   
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
    
    unless ($self->{$analysis.'_prediction'}) {
        my ($prediction, $score) = $self->_protein_function_prediction($analysis);
        $self->{$analysis}->{score} = $score;
        $self->{$analysis}->{prediction} = $prediction;
    }
    
    return $self->{$analysis}->{_prediction};
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

=head2 condel_prediction

  Description: Return the Condel (Consensus Deleteriousness) prediction for this allele that integrates
               the SIFT and Polyphen-2 scores
  Returntype : string (one of 'neutral', 'deleterious', 'non_computable_was') if this is a non-synonymous 
               mutation and predictions for this substitution are available from both SIFT and PolyPhen, 
               undef otherwise
  Exceptions : none
  Status     : At Risk

=cut

sub condel_prediction {
    my ($self, $condel_prediction) = @_;

    $self->{condel_prediction} = $condel_prediction if $condel_prediction;

    unless ($self->{condel_prediction}) {

        my $sift_score = $self->sift_score;
        my $pph_score  = $self->polyphen_score;
        my $pph_pred   = $self->polyphen_prediction;

        # we can only run Condel when we have predictions from both sift and polyphen

        if (defined $pph_score && defined $sift_score && ($pph_pred ne 'unknown') ) {
            my ($prediction, $score) = get_condel_prediction($sift_score, $pph_score);

            $self->{condel_prediction}  = $prediction;
            $self->{condel_score}       = $score;
        }
    }

    return $self->{condel_prediction};
}

=head2 condel_score

  Description: Return the Condel (Consensus Deleteriousness) score for this allele that integrates
               the SIFT and Polyphen-2 scores
  Returntype : float between 0 and 1 if this is a missense mutation and a prediction is 
               computable, -1 if SIFT and PolyPhen scores are available but Condel is unable 
               to compute a weighted average score, and undef otherwise
  Exceptions : none
  Status     : At Risk

=cut

sub condel_score {
    my ($self, $condel_score) = @_;

    $self->{condel_score} = $condel_score if defined $condel_score;

    # call condel_prediction to set the score if we need it
    $self->condel_prediction unless $self->{condel_score};

    return $self->{condel_score};
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
    return _hgvs_generic(@_,'coding');
}

=head2 hgvs_protein

  Description: Return a string representing the protein-level effect of this allele in HGVS format
  Returntype : string or undef if this allele is not in the CDS 
  Exceptions : none
  Status     : At Risk

=cut

sub hgvs_protein {
    return _hgvs_generic(@_,'protein');
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

1;
