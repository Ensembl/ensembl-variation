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
    
    unless ($self->{amino_acid}) {
        # just translate the codon sequence to establish the amino acid
        
        my $codon_seq = Bio::Seq->new(
            -seq        => $self->codon,
            -moltype    => 'dna',
            -alphabet   => 'dna',
        );
    
        $self->{amino_acid} = $codon_seq->translate->seq;
    }
    
    return $self->{amino_acid};
}

sub codon {
    my ($self, $codon) = @_;
    
    $self->{codon} = $codon if $codon;
    
    unless ($self->{codon}) {
        
        # calculate the codon sequence
    
        my $seq = $self->seq;
        
        if ($seq eq '-') {
            $self->{amino_acid} = $seq;
            $self->{codon} = $seq;
        }
        else {
            my $tv = $self->transcript_variation;         
            my $cds = $tv->translateable_seq;
            
            # calculate necessary coords and lengths
            
            my $codon_cds_start = $tv->pep_start * 3 - 2;
            my $codon_cds_end   = $tv->pep_end * 3;
            my $codon_len       = $codon_cds_end - $codon_cds_start + 1;
            my $vf_nt_len       = $tv->cds_end - $tv->cds_start + 1;
            my $allele_len      = length($seq);
            
            # splice the allele sequence into the CDS
        
            substr($cds, $tv->cds_start-1, $allele_len) = $seq;
            
            # and extract the codon sequence
            
            $self->{codon} = substr($cds, $codon_cds_start-1, $codon_len + ($allele_len - $vf_nt_len));
        }
    }
    
    return $self->{codon};
}

sub codon_position {
    my ($self, $codon_pos) = @_;
    
    $self->{codon_position} = $codon_pos if defined $codon_pos;
    
    unless ($self->{codon_position}) {
        # TODO: calculate the codon position (or store it in the DB?)
    }
    
    return $self->{codon_position};
}

sub affects_cds {
    my ($self) = @_;

    unless (defined $self->{affects_cds}) {
        
        my $tv = $self->transcript_variation;
       
        my @pep_coords = @{ $tv->pep_coords };
       
        if (@pep_coords != 1) {
            $self->{affects_cds} = 0;
        }
        elsif ($pep_coords[0]->isa('Bio::EnsEMBL::Mapper::Gap')) {
            $self->{affects_cds} = 0;
        }
        else {
            $self->{affects_cds} = 1;
        }
    }

    return $self->{affects_cds};
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
    
    $self->{qq{hgvs_$reference}} = $notation if defined $notation;
    
    unless ($self->{qq{hgvs_$reference}}) {
        # TODO: calculate the HGVS notation on-the-fly
    }
    
    return $self->{qq{hgvs_$reference}};
}

1;
