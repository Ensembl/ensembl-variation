package Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele;

use strict;
use warnings;

sub new {
    my ($class) = @_;
    return bless {}, $class;
}

sub new_fast {
    my ($class, $hashref) = @_;
    return bless $hashref, $class;
}

sub variation_feature_overlap {
    my ($self, $variation_feature_overlap) = @_;
    $self->{variation_feature_overlap} = $variation_feature_overlap if $variation_feature_overlap;
    return $self->{variation_feature_overlap};
}

sub seq {
    my ($self, $seq) = @_;
    $self->{seq} = $seq if $seq;
    return $self->{seq};
}

sub aa {
    my ($self, $aa) = @_;
    $self->{aa} = $aa if $aa;
    return $self->{aa};
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
        if (@pep_coords != 1) {
            $self->{affects_cds} = 0;
        }
        elsif ($pep_coords[0]->isa('Bio::EnsEMBL::Mapper::Gap')) {
            $self->{affects_cds} = 0;
        }
        else {

            my $tran = $vfo->feature;

            # by now we have established that the vf falls entirely in the CDS, 
            # so store the various start and end coords we need

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
            
            my $vf_nt_len = $vfo->cds_end - $vfo->cds_start + 1;
            
            # now work out the effect on the CDS and peptide sequence

            my $peptide = $tran->translate->seq;
            my $cds     = $tran->translateable_seq;

            # find the reference amino acid and codon

            my $var_pep_len = $vfo->pep_end - $vfo->pep_start + 1;

            my $ref_aa = substr($peptide, $vfo->pep_start - 1, $var_pep_len);
        
            $ref_aa = '-' unless $ref_aa;

            print "PEPTIDE:\n$peptide\nstart: ".$vfo->pep_start." end: ".$vfo->pep_end." len: ".$var_pep_len."\nREFERENCE AA: $ref_aa\n";

            $vfo->reference_allele->aa($ref_aa);

            my $codon_cds_start = $vfo->pep_start * 3 - 2;
            my $codon_cds_end   = $vfo->pep_end * 3;
            my $codon_len       = $codon_cds_end - $codon_cds_start + 1;

            my $ref_codon = substr($cds, $codon_cds_start-1, $codon_len);
            
            $ref_codon = '-' unless $ref_codon;
            
            $vfo->reference_allele->codon($ref_codon);
            
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
                
                $self->aa($new_aa);
            }
            else {
                $self->codon('-');
                $self->aa('-');
            }
          
            $self->{affects_cds} = 1;
        }
    }

    return $self->{affects_cds};
}

sub consequences {
    my ($self, $consequence) = @_;

    if ($consequence) {
        my $consequences = $self->{consequences} ||= [];
        push @$consequences, $consequence; 
    }
    
    return $self->{consequences};
}

sub is_reference {
    my ($self, $is_reference) = @_;
    $self->{is_reference} = $is_reference if defined $is_reference;
    return $self->{is_reference};
}

sub calc_consequences {
    my ($self, $consequences) = @_;

    for my $consequence (@$consequences) {
        if ($consequence->predicate->($self)) {
            $self->consequences($consequence);
        }
    }
}

1;
