=head1 LICENSE

 Copyright (c) 1999-2011 The European Bioinformatics Institute and
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

package Bio::EnsEMBL::Variation::TranscriptVariationAllele;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele);

sub new_fast {
    my ($self, $hashref) = @_;
    
    # swap a transcript_variation argument for a variation_feature_overlap one

    if ($hashref->{transcript_variation}) {
        $hashref->{variation_feature_overlap} = delete $hashref->{transcript_variation};
    }
    
    # and call the superclass

    return $self->SUPER::new_fast($hashref);
}

sub transcript_variation {
    my $self = shift;
    return $self->variation_feature_overlap(@_);
}

sub transcript {
    my $self = shift;
    return $self->transcript_variation->transcript;
}

sub variation_feature {
    my $self = shift;
    return $self->transcript_variation->variation_feature;
}

sub pep_allele_string {
    my ($self) = @_;
    
    my $pep = $self->peptide;
    
    return undef unless $pep;
    
    my $ref_pep = $self->transcript_variation->get_reference_TranscriptVariationAllele->peptide;
    
    return $ref_pep ne $pep ? $ref_pep.'/'.$pep : $pep;
}

sub codon_allele_string {
    my ($self) = @_;
    
    my $codon = $self->codon;
    
    return undef unless $codon;
    
    my $ref_codon = $self->transcript_variation->get_reference_TranscriptVariationAllele->codon;
    
    return $ref_codon.'/'.$codon;
}

sub peptide {
    my ($self, $peptide) = @_;
    
    $self->{peptide} = $peptide if $peptide;
    
    unless ($self->{peptide}) {
        
        if (my $codon = $self->codon) {
            
            # the codon method can set the peptide in some circumstances 
            # so check here before we try an (expensive) translation
            return $self->{peptide} if $self->{peptide};
            
            # translate the codon sequence to establish the peptide allele
            
            # for mithocondrial dna we need to to use a different codon table
            my $codon_table = $self->transcript_variation->_codon_table;
            
            my $codon_seq = Bio::Seq->new(
                -seq        => $codon,
                -moltype    => 'dna',
                -alphabet   => 'dna',
            );
        
            my $pep = $codon_seq->translate(undef, undef, undef, $codon_table)->seq;
            
            if (length($pep) < 1) {
                if (length($codon) % 3) {
                    # partial codon
                    $pep = 'X';
                }
                else {
                    $pep = '-';
                }
            }
           
            $self->{peptide} = $pep;
        }
    }
    
    return $self->{peptide};
}

sub codon {
    my ($self, $codon) = @_;
    
    $self->{codon} = $codon if defined $codon;
    
    my $tv = $self->transcript_variation;      
    
    return undef unless $tv->translation_start;
    
    return undef if $self->variation_feature_seq =~ /[^ACGT\-]/i;
    
    unless ($self->{codon}) {
      
        # try to calculate the codon sequence
    
        my $seq = $self->feature_seq;
        
        $seq = '' if $seq eq '-';
        
        # calculate necessary coords and lengths
        
        my $codon_cds_start = $tv->translation_start * 3 - 2;
        my $codon_cds_end   = $tv->translation_end * 3;
        my $codon_len       = $codon_cds_end - $codon_cds_start + 1;
        my $vf_nt_len       = $tv->cds_end - $tv->cds_start + 1;
        my $allele_len      = length($seq);
        
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

sub polyphen_prediction {
    my ($self, $polyphen_prediction) = @_;
    
    $self->{polyphen_prediction} = $polyphen_prediction if $polyphen_prediction;
    
    $self->{polyphen_prediction} = $self->_nsSNP_prediction('polyphen')
        unless $self->{polyphen_prediction};
    
    return $self->{polyphen_prediction};
}

sub sift_prediction {
    my ($self, $sift_prediction) = @_;
    
    $self->{sift_prediction} = $sift_prediction if $sift_prediction;
    
    $self->{sift_prediction} = $self->_nsSNP_prediction('sift')
        unless $self->{sift_prediction};
    
    return $self->{sift_prediction};
}

sub _nsSNP_prediction {
    my ($self, $program) = @_;

    # we can only get results for variants that cause a single amino acid substitution, 
    # so check the peptide allele string first

    if ($self->pep_allele_string && $self->pep_allele_string =~ /^[A-Z]\/[A-Z]$/) {
        if (my $adap = $self->transcript_variation->{adaptor}) {
            return $adap->_get_nsSNP_prediction($program, $self);
        }
    }
    
    return undef;
}

sub hgvs_genomic {
    return _hgvs_generic(@_,'genomic');
}
sub hgvs_coding {
    return _hgvs_generic(@_,'coding');
}
sub hgvs_protein {
    return _hgvs_generic(@_,'protein');
}
sub hgvs_rna {
    return _hgvs_generic(@_,'rna');
}
sub hgvs_mitochondrial {
    return _hgvs_generic(@_,'mitochondrial');
}

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
