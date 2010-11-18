package Bio::EnsEMBL::Variation::TranscriptVariationNew;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::VariationFeatureOverlap);

#use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(overlap);

#use Inline C => <<'END_C';
#
#int overlap (int f1_start, int f1_end, int f2_start, int f2_end) {
#    return (f1_end >= f2_start && f1_start <= f2_end);
#}
#
#END_C

sub overlap {
    my ( $f1_start, $f1_end, $f2_start, $f2_end ) = @_;
    return ($f1_end >= $f2_start and $f1_start <= $f2_end);
}

sub transcript {
    my $self = shift;
    
    if (my $tran_id = $self->{_feature_stable_id}) {
        
        # lazy-load the Transcript
        
        if (my $adap = $self->{adaptor}) {
            if (my $ta = $adap->db->dnadb->get_TranscriptAdaptor) {
                if (my $tran = $ta->fetch_by_stable_id($tran_id)) {
                    $self->{feature} = $tran;
                    delete $self->{_feature_stable_id};
                }
            }
        }
    }
    
    return $self->feature(@_);
}

sub cdna_start {
    my ($self, $cdna_start) = @_;
    
    $self->{cdna_start} = $cdna_start if $cdna_start;
    
    unless ($self->{cdna_start}) {
        my $cdna_coords = $self->cdna_coords;
        
        return undef if @$cdna_coords != 1;
        
        $self->{cdna_start} = $cdna_coords->[0]->start;
        $self->{cdna_end} = $cdna_coords->[0]->end;
    }
    
    return $self->{cdna_start};
}

sub cdna_end {
    my ($self, $cdna_end) = @_;
    
    $self->{cdna_end} = $cdna_end if $cdna_end;
    
    $self->cdna_start unless $self->{cdna_end};
    
    return $self->{cdna_end};
}

sub cds_start {
    my ($self, $cds_start) = @_;
    
    $self->{cds_start} = $cds_start if $cds_start;
    
    unless ($self->{cds_start}) {
        my $cds_coords = $self->cds_coords;
        
        return undef if @$cds_coords != 1;
        
        my $exon_phase = $self->transcript->start_Exon->phase;
        
        $self->{cds_start} = $cds_coords->[0]->start + ($exon_phase > 0 ? $exon_phase : 0);
        $self->{cds_end} = $cds_coords->[0]->end + ($exon_phase > 0 ? $exon_phase : 0);;
    }
    
    return $self->{cds_start};
}

sub cds_end {
    my ($self, $cds_end) = @_;
    
    $self->{cds_end} = $cds_end if $cds_end;
    
    $self->cds_start unless $self->{cds_end};
    
    return $self->{cds_end};
}

sub pep_start {
    my ($self, $pep_start) = @_;
    
    $self->{pep_start} = $pep_start if $pep_start;
    
    unless ($self->{pep_start}) {
        my $pep_coords = $self->pep_coords;
        
        return undef if @$pep_coords != 1;
        
        $self->{pep_start} = $pep_coords->[0]->start;
        $self->{pep_end} = $pep_coords->[0]->end;
    }
    
    return $self->{pep_start};
}

sub pep_end {
    my ($self, $pep_end) = @_;
    
    $self->{pep_end} = $pep_end if $pep_end;
    
    $self->pep_start unless $self->{pep_end};
    
    return $self->{pep_end};
}

sub cdna_coords {
    my ($self) = @_;
    
    unless ($self->{cdna_coords}) {
        my $vf   = $self->variation_feature;
        my $tran = $self->feature; 
        $self->{cdna_coords} = [ $self->mapper->genomic2cdna($vf->start, $vf->end, $tran->strand) ];
    }
    
    return $self->{cdna_coords};
}

sub cds_coords {
    my ($self) = @_;
    
    unless ($self->{cds_coords}) {
        my $vf   = $self->variation_feature;
        my $tran = $self->feature; 
        $self->{cds_coords} = [ $self->mapper->genomic2cds($vf->start, $vf->end, $tran->strand) ];
    }
    
    return $self->{cds_coords};
}

sub pep_coords {
    my ($self) = @_;
    
    unless ($self->{pep_coords}) {
        my $vf   = $self->variation_feature;
        my $tran = $self->feature; 
        $self->{pep_coords} = [ $self->mapper->genomic2pep($vf->start, $vf->end, $tran->strand) ];
    }
    
    return $self->{pep_coords};
}

sub intron_effects {
    my ($self) = @_;
    
    # this method is a major bottle neck in the effect calculation code so 
    # we cache results and use local variables instead of method calls where
    # possible to speed things up - caveat bug-fixer!
    
    unless ($self->{intron_effects}) {
        
        my $vf = $self->variation_feature;
        
        my $intron_effects = {};
        
        my $found_effect = 0;
        
        my $vf_start = $vf->start;
        my $vf_end   = $vf->end;

        for my $intron (@{ $self->introns }) {
            
            my $intron_start = $intron->start;
            my $intron_end   = $intron->end;
            
            # under various circumstances the genebuild process can introduce
            # artificial short (<= 12 nucleotide) introns into transcripts
            # (e.g. to deal with errors in the reference sequence etc.), we
            # don't want to categorise variations that fall in these introns
            # as intronic, or as any kind of splice variant
            
            my $frameshift_intron = ( abs($intron_end - $intron_start) <= 12 );
            
            # the order of these checks is deliberately designed to minimise the number
            # of calls we make to overlap because we can sometimes establish several
            # intron effects within one call and then break out of the loop with last
            
            if (overlap($vf_start, $vf_end, $intron_start, $intron_start+2)) {
                
                if ($frameshift_intron) {
                    $intron_effects->{within_frameshift_intron} = 1;
                }
                else {
                    $intron_effects->{start_splice_site} = 1;
                    $intron_effects->{splice_region} = 1;
                    $intron_effects->{intronic} = 1;
                }
                
                last;
            }
            
            if (overlap($vf_start, $vf_end, $intron_end-2, $intron_end)) {
                
                if ($frameshift_intron) {
                    $intron_effects->{within_frameshift_intron} = 1;
                }
                else {
                    $intron_effects->{end_splice_site} = 1;
                    $intron_effects->{splice_region} = 1;
                    $intron_effects->{intronic} = 1;
                }
                
                last;
            }
            
            # we don't last out of this check, because a variation can be both
            # intronic and splice_region, so we just set a flag
            
            if (overlap($vf_start, $vf_end, $intron_start, $intron_end)) {
                
                if ($frameshift_intron) {
                    $intron_effects->{within_frameshift_intron} = 1;
                }
                else {
                    $intron_effects->{intronic} = 1;
                }
                
                $found_effect = 1;
            }
            
            if ( ( overlap($vf_start, $vf_end, $intron_start-3, $intron_start+8) or
                   overlap($vf_start, $vf_end, $intron_end-8, $intron_end+3) ) and 
                   not $frameshift_intron ) {
                   
                $intron_effects->{splice_region} = 1;
                
                $found_effect = 1;
            }
                
            last if $found_effect;
        }
        
        $self->{intron_effects} = $intron_effects;       
    }

    return $self->{intron_effects};
}

# NB: the methods below all cache their data in the associated transcript itself, this
# gives a significant speed up when you are calculating the effect of all variations
# on a transcript, and means that the cache will be freed when the transcript itself
# is garbage collected rather than us having to maintain a transcript feature cache 
# ourselves

sub introns {
    my ($self) = @_;
    
    my $tran = $self->feature;
    
    my $introns = $tran->{_variation_effect_feature_cache}->{introns} ||= $tran->get_all_Introns;
    
    return $introns;
}

sub translateable_seq {
    my ($self) = @_;
    
    my $tran = $self->feature;
    
    my $tran_seq = $tran->{_variation_effect_feature_cache}->{translateable_seq} ||= $tran->translateable_seq;
    
    return $tran_seq;
}

sub mapper {
    my ($self) = @_;
    
    my $tran = $self->feature;
    
    my $mapper = $tran->{_variation_effect_feature_cache}->{mapper} ||= $tran->get_TranscriptMapper;
    
    return $mapper;
}
  
sub peptide {
    my ($self) = @_;
    
    my $tran = $self->feature;
    
    my $peptide = $tran->{_variation_effect_feature_cache}->{peptide};
    
    unless ($peptide) {
        my $translation = $tran->translate;
        $peptide = $translation ? $translation->seq : undef;
        $tran->{_variation_effect_feature_cache}->{peptide} = $peptide;
    }
    
    return $peptide;
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
    my $hgvs = shift;
    
    #ÊThe rna and mitochondrial modes have not yet been implemented, so return undef in case we get a call to these
    return undef if ($reference =~ m/rna|mitochondrial/);
    
    # The HGVS subroutine
    my $sub = qq{hgvs_$reference};
    
    # Loop over the TranscriptVariationAllele objects associated with this TranscriptVariation
    foreach my $tv_allele (@{$self->alt_alleles()}) {
        
        #ÊIf an HGVS hash was supplied and the allele exists as key, set the HGVS notation for this allele
        if (defined($hgvs) && (my $notation = $hgvs->{$tv_allele->seq()})) {
            $tv_allele->$sub($notation);
        }
        # Else, add the HGVS notation for this allele to the HGVS hash
        else {
            $hgvs->{$tv_allele->seq()} = $tv_allele->$sub();
        }
    }
    
    return $hgvs;
}

1;
