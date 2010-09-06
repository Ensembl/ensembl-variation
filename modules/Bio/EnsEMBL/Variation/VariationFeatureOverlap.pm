package Bio::EnsEMBL::Variation::VariationFeatureOverlap;

use strict;
use warnings;

#use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(overlap);

use Inline C => <<'END_C';

int overlap (int f1_start, int f1_end, int f2_start, int f2_end) {
    return (f1_end >= f2_start && f1_start <= f2_end);
}

END_C

sub new {
    my ($class) = @_;
    return bless {}, $class;
}

sub new_fast {
    my ($class, $hashref) = @_;
    return bless $hashref, $class;
}

sub variation_feature {
    my ($self, $variation_feature) = @_;
    $self->{variation_feature} = $variation_feature if $variation_feature;
    return $self->{variation_feature};
}

sub feature {
    my ($self, $feature) = @_;
    $self->{feature} = $feature if $feature;
    return $self->{feature};
}

sub feature_type {
    my ($self, $feature_type) = @_;
    $self->{feature_type} = $feature_type if $feature_type;
    return $self->{feature_type};
}

sub cdna_start {
    my ($self, $cdna_start) = @_;
    $self->{cdna_start} = $cdna_start if $cdna_start;
    return $self->{cdna_start};
}

sub cdna_end {
    my ($self, $cdna_end) = @_;
    $self->{cdna_end} = $cdna_end if $cdna_end;
    return $self->{cdna_end};
}

sub cds_start {
    my ($self, $cds_start) = @_;
    $self->{cds_start} = $cds_start if $cds_start;
    return $self->{cds_start};
}

sub cds_end {
    my ($self, $cds_end) = @_;
    $self->{cds_end} = $cds_end if $cds_end;
    return $self->{cds_end};
}

sub pep_start {
    my ($self, $pep_start) = @_;
    $self->{pep_start} = $pep_start if $pep_start;
    return $self->{pep_start};
}

sub pep_end {
    my ($self, $pep_end) = @_;
    $self->{pep_end} = $pep_end if $pep_end;
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

sub reference_allele {
    my ($self, $reference_allele) = @_;
    $self->{reference_allele} = $reference_allele if $reference_allele;
    return $self->{reference_allele};
}

sub alt_alleles {
    my ($self, @new_alt_alleles) = @_;
    
    $self->{alt_alleles} ||= [];

    if (@new_alt_alleles) {
        my $alt_alleles = $self->{alt_alleles};
        push @$alt_alleles, @new_alt_alleles;
    }

    return $self->{alt_alleles};
}

sub alleles {
    my ($self) = @_;

    my $alleles = $self->alt_alleles;
    
    unshift @$alleles, $self->reference_allele;

    return $alleles;
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

sub introns {
    my ($self, $introns) = @_;
    $self->{introns} = $introns if $introns;
    return $self->{introns};
}

sub translateable_seq {
    my ($self, $translateable_seq) = @_;
    $self->{translateable_seq} = $translateable_seq if $translateable_seq;
    return $self->{translateable_seq};
}

sub peptide {
    my ($self, $peptide) = @_;
    $self->{peptide} = $peptide if $peptide;
    return $self->{peptide};
}

sub mapper {
    my ($self, $mapper) = @_;
    $self->{mapper} = $mapper if $mapper;
    return $self->{mapper};
}

1;


