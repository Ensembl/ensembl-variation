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
        # translate the codon sequence to establish the amino acid
        
        # for mithocondrial dna we need to to use a different codon table
        my $codon_table = $self->transcript_variation->codon_table;
        
        my $codon_seq = Bio::Seq->new(
            -seq        => $self->codon,
            -moltype    => 'dna',
            -alphabet   => 'dna',
        );
    
        $self->{amino_acid} = $codon_seq->translate(undef, undef, undef, $codon_table)->seq;
    }
    
    return $self->{amino_acid};
}

sub codon {
    my ($self, $codon) = @_;
    
    $self->{codon} = $codon if $codon;
    
    unless ($self->{codon}) {
        
        # calculate the codon sequence
    
        my $seq = $self->feature_seq;
        
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

sub consequences {
    my ($self, @consequences) = @_;
    
    my $cons = $self->{consequences} ||= [];
    
    if (@consequences) {
        push @$cons, @consequences;
    }
    
    unless (@$cons) {
        
        # calculate consequences on the fly
        
        if (my $overlap_cons = $self->transcript_variation->overlap_consequences) {
            for my $oc (@$overlap_cons) {
                if ($oc->predicate->($self)) {
                    push @$cons, $oc;
                }
            }
        }
    }
    
    return $cons;
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
    
    #ÊThe rna and mitochondrial modes have not yet been implemented, so return undef in case we get a call to these
    return undef if ($reference =~ m/rna|mitochondrial/);
    
    my $sub = qq{hgvs_$reference};
    
    $self->{$sub} = $notation if defined $notation;
    
    unless ($self->{$sub}) {
        # Calculate the HGVS notation on-the-fly and pass it to the TranscriptVariation in order to distribute the result to the other alleles
        my $reference_feature;
        $reference_feature = $self->transcript unless ($reference eq 'genomic');
        $self->transcript_variation->$sub($self->variation_feature->get_all_hgvs_notations($reference_feature,substr($reference,0,1)));
    }
    
    return $self->{$sub};
}

1;
