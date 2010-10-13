package Bio::EnsEMBL::Variation::TranscriptVariationAllele;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele);

sub transcript_variation {
    my $self = shift;
    return $self->variation_feature_overlap;
}

sub transcript {
    my $self = shift;
    return $self->variation_feature_overlap->feature;
}

sub variation_feature {
    my $self = shift;
    return $self->variation_feature_overlap->variation_feature;
}

sub amino_acid {
    my ($self, $amino_acid) = @_;
    $self->{amino_acid} = $amino_acid if $amino_acid;
    return $self->{amino_acid};
}

sub codon {
    my ($self, $codon) = @_;
    $self->{codon} = $codon if $codon;
    return $self->{codon};
}

sub affects_cds {
    my ($self) = @_;

    unless (defined $self->{affects_cds}) {
        
        my $vfo = $self->variation_feature_overlap;
       
        my @pep_coords = @{ $vfo->pep_coords };
       
        if (not $vfo->feature->isa('Bio::EnsEMBL::Transcript')) {
            $self->{affects_cds} = 0;
        }
        elsif (@pep_coords != 1) {
            $self->{affects_cds} = 0;
        }
        elsif ($pep_coords[0]->isa('Bio::EnsEMBL::Mapper::Gap')) {
            $self->{affects_cds} = 0;
        }
        else {

            # by now we have established that the vf falls entirely in the CDS, 
            # so store the various start and end coords we need

            my $tran = $vfo->feature;
            
            my ($cds, $codon_cds_start, $codon_cds_end, $codon_len);
            
            # there's some stuff we only need to calculate once per vf, rather
            # than once per allele
            
            if (defined $vfo->pep_start) {
                
                # but we do need these in any case
                
                $codon_cds_start = $vfo->pep_start * 3 - 2;
                $codon_cds_end   = $vfo->pep_end * 3;
                $codon_len       = $codon_cds_end - $codon_cds_start + 1;
                $cds             = $vfo->translateable_seq;
            }
            else {
                
                # this is the stuff we only need calculate once per vf
                
                my $exon_phase = $tran->start_Exon->phase;
    
                my $pep_coord = $vfo->pep_coords->[0];
                $vfo->pep_start($pep_coord->start);
                $vfo->pep_end($pep_coord->end);
    
                my $cds_coord = $vfo->cds_coords->[0];
                $vfo->cds_start($cds_coord->start + ($exon_phase > 0 ? $exon_phase : 0));
                $vfo->cds_end($cds_coord->end + ($exon_phase > 0 ? $exon_phase : 0));
    
                my $cdna_coord = $vfo->cdna_coords->[0];
                $vfo->cdna_start($cdna_coord->start);
                $vfo->cdna_end($cdna_coord->end);
                
                # find the reference amino acid and codon
    
                my $peptide = $vfo->peptide;
                $cds        = $vfo->translateable_seq;
 
                my $var_pep_len = $vfo->pep_end - $vfo->pep_start + 1;
    
                my $ref_aa = substr($peptide, $vfo->pep_start - 1, $var_pep_len);
            
                $ref_aa = '-' unless $ref_aa;
    
                #print "PEPTIDE:\n$peptide\nstart: ".$vfo->pep_start." end: ".$vfo->pep_end." len: ".$var_pep_len."\nREFERENCE AA: $ref_aa\n";
    
                $vfo->reference_allele->amino_acid($ref_aa);
    
                $codon_cds_start = $vfo->pep_start * 3 - 2;
                $codon_cds_end   = $vfo->pep_end * 3;
                $codon_len       = $codon_cds_end - $codon_cds_start + 1;
    
                my $ref_codon = substr($cds, $codon_cds_start-1, $codon_len);
                
                $ref_codon = '-' unless $ref_codon;
                
                $vfo->reference_allele->codon($ref_codon);
            }
            
            # everything else does need to be calculated for each allele
            
            my $vf_nt_len = $vfo->cds_end - $vfo->cds_start + 1;
            
            # calculate the effect this allele has on the reference

            my $seq = $self->seq;

            $seq =~ s/-//;
            
            if ($seq) {
                my $new_cds = $cds;
                
                my $allele_len = length($seq);
    
                substr($new_cds, $vfo->cds_start-1, $allele_len) = $seq;
                
                my $new_codon = substr($new_cds, $codon_cds_start-1, $codon_len + ($allele_len - $vf_nt_len));
    
                $self->codon($new_codon);
    
                my $new_codon_seq = Bio::Seq->new(
                    -seq        => $new_codon,
                    -moltype    => 'dna',
                    -alphabet   => 'dna'
                );
    
                my $new_aa = $new_codon_seq->translate->seq;
    
                $new_aa = '-' unless $new_aa;
                
                $self->amino_acid($new_aa);
            }
            else {
                $self->codon('-');
                $self->amino_acid('-');
            }
          
            $self->{affects_cds} = 1;
        }
    }

    return $self->{affects_cds};
}

1;
