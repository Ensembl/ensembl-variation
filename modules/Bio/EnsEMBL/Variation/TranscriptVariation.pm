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

package Bio::EnsEMBL::Variation::TranscriptVariation;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::TranscriptVariationAllele;
use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(overlap within_cds);

use base qw(Bio::EnsEMBL::Variation::VariationFeatureOverlap);

sub new {
    my $class = shift;

    my %args = @_;

    # swap a '-transcript' argument for a '-feature' one for the superclass

    for my $arg (keys %args) {
        if (lc($arg) eq '-transcript') {
            $args{'-feature'} = delete $args{$arg};
        }
    }

    # call the superclass constructor
    my $self = $class->SUPER::new(%args) || return undef;

    # rebless the alleles from vfoas to tvas
    map { bless $_, 'Bio::EnsEMBL::Variation::TranscriptVariationAllele' } @{ $self->alleles };
    
    return $self;
}

sub transcript_stable_id {
    my $self = shift;
    return $self->SUPER::feature_stable_id(@_);
}

sub transcript {
    my ($self, $transcript) = @_;
    return $self->SUPER::feature($transcript, 'Transcript');
}

sub codon_position {
    my ($self, $codon_pos) = @_;
    
    $self->{codon_position} = $codon_pos if defined $codon_pos;
    
    unless ($self->{codon_position}) {
        $self->{codon_position} = (($self->cdna_start - $self->transcript->cdna_coding_start) % 3) + 1;
    }
    
    return $self->{codon_position};
}

sub cdna_start {
    my ($self, $cdna_start) = @_;
    
    $self->{cdna_start} = $cdna_start if $cdna_start;
    
    unless ($self->{cdna_start}) {
        my $cdna_coords = $self->cdna_coords;
        
        return undef if (@$cdna_coords != 1 || $cdna_coords->[0]->isa('Bio::EnsEMBL::Mapper::Gap'));
        
        $self->{cdna_start} = $cdna_coords->[0]->start;
        $self->{cdna_end}   = $cdna_coords->[0]->end;
    }
    
    return $self->{cdna_start};
}

sub cdna_end {
    my ($self, $cdna_end) = @_;
    
    $self->{cdna_end} = $cdna_end if $cdna_end;
    
    # call cdna_start to calculate the start and end
    $self->cdna_start unless $self->{cdna_end};
    
    return $self->{cdna_end};
}

sub cds_start {
    my ($self, $cds_start) = @_;
    
    $self->{cds_start} = $cds_start if $cds_start;
    
    unless ($self->{cds_start}) {
        my $cds_coords = $self->cds_coords;
        
        return undef if (@$cds_coords != 1 || $cds_coords->[0]->isa('Bio::EnsEMBL::Mapper::Gap'));
        
        my $exon_phase = $self->transcript->start_Exon->phase;
        
        $self->{cds_start} = $cds_coords->[0]->start + ($exon_phase > 0 ? $exon_phase : 0);
        $self->{cds_end}   = $cds_coords->[0]->end   + ($exon_phase > 0 ? $exon_phase : 0);
    }
    
    return $self->{cds_start};
}

sub cds_end {
    my ($self, $cds_end) = @_;
    
    $self->{cds_end} = $cds_end if $cds_end;
    
    # call cds_start to calculate the start and end
    $self->cds_start unless $self->{cds_end};
    
    return $self->{cds_end};
}

sub translation_start {
    my ($self, $translation_start) = @_;
    
    $self->{translation_start} = $translation_start if $translation_start;
    
    unless ($self->{translation_start}) {
        my $translation_coords = $self->translation_coords;
        
        return undef if (@$translation_coords != 1 || $translation_coords->[0]->isa('Bio::EnsEMBL::Mapper::Gap'));
        
        $self->{translation_start} = $translation_coords->[0]->start;
        $self->{translation_end}   = $translation_coords->[0]->end;
    }
    
    return $self->{translation_start};
}

sub translation_end {
    my ($self, $translation_end) = @_;
    
    $self->{translation_end} = $translation_end if $translation_end;
    
    # call translation_start to calculate the start and end
    $self->translation_start unless $self->{translation_end};
    
    return $self->{translation_end};
}

sub cdna_coords {
    my ($self) = @_;
    
    unless ($self->{cdna_coords}) {
        my $vf   = $self->variation_feature;
        my $tran = $self->transcript; 
        $self->{cdna_coords} = [ $self->mapper->genomic2cdna($vf->seq_region_start, $vf->seq_region_end, $tran->strand) ];
    }
    
    return $self->{cdna_coords};
}

sub cds_coords {
    my ($self) = @_;
    
    unless ($self->{cds_coords}) {
        my $vf   = $self->variation_feature;
        my $tran = $self->transcript; 
        $self->{cds_coords} = [ $self->mapper->genomic2cds($vf->seq_region_start, $vf->seq_region_end, $tran->strand) ];
    }
    
    return $self->{cds_coords};
}

sub translation_coords {
    my ($self) = @_;
    
    unless ($self->{translation_coords}) {
        my $vf   = $self->variation_feature;
        my $tran = $self->transcript; 
        $self->{translation_coords} = [ $self->mapper->genomic2pep($vf->seq_region_start, $vf->seq_region_end, $tran->strand) ];
    }
    
    return $self->{translation_coords};
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
        
        my $vf_start = $vf->seq_region_start;
        my $vf_end   = $vf->seq_region_end;

        for my $intron (@{ $self->introns }) {
            
            my $intron_start = $intron->seq_region_start;
            my $intron_end   = $intron->seq_region_end;
            
            # under various circumstances the genebuild process can introduce
            # artificial short (<= 12 nucleotide) introns into transcripts
            # (e.g. to deal with errors in the reference sequence etc.), we
            # don't want to categorise variations that fall in these introns
            # as intronic, or as any kind of splice variant
            
            my $frameshift_intron = ( abs($intron_end - $intron_start) <= 12 );
            
            if ($frameshift_intron) {
                if (overlap($vf_start, $vf_end, $intron_start, $intron_end)) {
                    $intron_effects->{within_frameshift_intron} = 1;
                    last;
                }
            }
            
            # the order of these checks is deliberately designed to minimise the number
            # of calls we make to overlap because we can sometimes establish several
            # intron effects within one call and then break out of the loop with last
            
            if (overlap($vf_start, $vf_end, $intron_start, $intron_start+1)) {
                $intron_effects->{start_splice_site} = 1;
                #$intron_effects->{splice_region} = 1;
                #$intron_effects->{intronic} = 1;
                
                last;
            }
            
            if (overlap($vf_start, $vf_end, $intron_end-1, $intron_end)) {
                $intron_effects->{end_splice_site} = 1;
                #$intron_effects->{splice_region} = 1;
                #$intron_effects->{intronic} = 1;
                
                last;
            }
            
            # we don't last out of this check, because a variation can be both
            # intronic and splice_region, so we just set a flag
            
            if (overlap($vf_start, $vf_end, $intron_start, $intron_end)) {
                $intron_effects->{intronic} = 1;
                
                $found_effect = 1;
            }
            
            # the definition of splice_region (SO:0001630) is "within 1-3 bases 
            # of the exon or 3-8 bases of the intron", the intron start is the 
            # first base of the intron so we only need to add or subtract 7 from 
            # it to get the correct coordinate (we don't check if it falls in the
            # donor or acceptor splice sites because that would have been 
            # established above)
            
            if ( ( overlap($vf_start, $vf_end, $intron_start-3, $intron_start+7) or
                   overlap($vf_start, $vf_end, $intron_end-7,   $intron_end+3) ) and 
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
    
    my $tran = $self->transcript;
    
    my $introns = $tran->{_variation_effect_feature_cache}->{introns} ||= $tran->get_all_Introns;
    
    return $introns;
}

sub translateable_seq {
    my ($self) = @_;
    
    my $tran = $self->transcript;
    
    my $tran_seq = $tran->{_variation_effect_feature_cache}->{translateable_seq} ||= $tran->translateable_seq;
    
    return $tran_seq;
}

sub mapper {
    my ($self) = @_;
    
    my $tran = $self->transcript;
    
    my $mapper = $tran->{_variation_effect_feature_cache}->{mapper} ||= $tran->get_TranscriptMapper;
    
    return $mapper;
}
  
sub peptide {
    my ($self) = @_;
    
    my $tran = $self->transcript;
    
    my $peptide = $tran->{_variation_effect_feature_cache}->{peptide};
    
    unless ($peptide) {
        my $translation = $tran->translate;
        $peptide = $translation ? $translation->seq : undef;
        $tran->{_variation_effect_feature_cache}->{peptide} = $peptide;
    }
    
    return $peptide;
}

sub codon_table {
    my ($self) = @_;
    
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
        if (defined($hgvs)) {
            my $notation = $hgvs->{$tv_allele->variation_feature_seq()};
            $tv_allele->$sub($notation) if defined $notation;
        }
        # Else, add the HGVS notation for this allele to the HGVS hash
        else {
            $hgvs->{$tv_allele->variation_feature_seq()} = $tv_allele->$sub();
        }
    }
    
    return $hgvs;
}

sub pep_allele_string {
    my $self = shift;
    
    unless ($self->{_pep_allele_string}) {
        
        my @peptides = grep {defined $_} map { $_->peptide } @{ $self->alleles };
        
        $self->{_pep_allele_string} = join '/', @peptides;
    }
    
    return $self->{_pep_allele_string};
}

sub codons {
    my $self = shift;
    
    unless ($self->{_codon_allele_string}) {
        
        my @codons = grep {defined $_} map { $_->codon } @{ $self->alleles };
        
        $self->{_codon_allele_string} = join '/', @codons;
    }
    
    return $self->{_codon_allele_string};
}

sub affects_transcript {
    my $self = shift;
    return scalar grep { within_cds($_) } @{ $self->alt_alleles }; 
}

1;
