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

sub pep_allele_string {
    my ($self) = @_;
    
    my $pep = $self->peptide;
    
    return undef unless $pep;
    
    my $ref_pep = $self->transcript_variation->reference_allele->peptide;
    
    return $ref_pep ne $pep ? $ref_pep.'/'.$pep : $pep;
}

sub codon_allele_string {
    my ($self) = @_;
    
    my $codon = $self->codon;
    
    return undef unless $codon;
    
    my $ref_codon = $self->transcript_variation->reference_allele->codon;
    
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
            my $codon_table = $self->transcript_variation->codon_table;
            
            my $codon_seq = Bio::Seq->new(
                -seq        => $codon,
                -moltype    => 'dna',
                -alphabet   => 'dna',
            );
        
            my $pep = $codon_seq->translate(undef, undef, undef, $codon_table)->seq;
            
            $pep = '-' if length($pep) < 1;
            
            $self->{peptide} = $pep;
        }
    }
    
    return $self->{peptide};
}

sub codon {
    my ($self, $codon) = @_;
    
    $self->{codon} = $codon if $codon;
    
    my $tv = $self->transcript_variation;      
    
    return undef unless $tv->pep_start;
    
    return undef if $self->variation_feature_seq =~ /[^ACGT\-]/i;
    
    unless ($self->{codon}) {
      
        # calculate the codon sequence
    
        my $seq = $self->feature_seq;
        
        $seq = '' if $seq eq '-';
        
        my $cds = $tv->translateable_seq;
        
        # calculate necessary coords and lengths
        
        my $codon_cds_start = $tv->pep_start * 3 - 2;
        my $codon_cds_end   = $tv->pep_end * 3;
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

sub consequence_types {
    my ($self, @new_consequences) = @_;
    
    my $cons = $self->{consequence_types};
    
    if (@new_consequences) {
        $cons ||= [];
        push @$cons, @new_consequences;
    }
    
    unless (defined $cons) {
        
        print "OTF consequences...\n";
        
        # calculate consequences on the fly
        
        $cons = [];
        
        if (my $overlap_cons = $self->transcript_variation->overlap_consequences) {
            for my $oc (@$overlap_cons) {
                if ($oc->predicate->($self)) {
                    push @$cons, $oc;
                }
            }
        }
    }
    
    $self->{consequence_types} = $cons;
    
    return $cons;
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
