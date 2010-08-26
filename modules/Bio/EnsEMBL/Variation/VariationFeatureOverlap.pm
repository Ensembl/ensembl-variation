package Bio::EnsEMBL::Variation::VariationFeatureOverlap;

use strict;
use warnings;

sub new {
    my ($class) = @;
    return bless {}, $class;
}

sub new_fast {
    my ($class, $hashref) = @;
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
    my ($self, $cdna_coords) = @_;
    $self->{cdna_coords} = $cdna_coords if $cdna_coords;
    return $self->{cdna_coords};
}

sub cds_coords {
    my ($self, $cds_coords) = @_;
    $self->{cds_coords} = $cds_coords if $cds_coords;
    return $self->{cds_coords};
}

sub pep_coords {
    my ($self, $pep_coords) = @_;
    $self->{pep_coords} = $pep_coords if $pep_coords;
    return $self->{pep_coords};
}

sub reference_allele {
    my ($self, $reference_allele) = @_;
    $self->{reference_allele} = $reference_allele if $reference_allele;
    return $self->{reference_allele};
}

sub alleles {
    my ($self, $allele) = @_;

    if ($allele) {
        my $alleles = $self->{alleles} ||= [];
        push @$alleles, $allele
    }

    return $self->{alleles};
}

1;

package Bio::EnsEMBL::Variation::OverlapConsequence;

use strict;
use warnings;

sub new {
    my ($class) = @;
    return bless {}, $class;
}

sub new_fast {
    my ($class, $hashref) = @;
    return bless $hashref, $class;
}

sub SO_term {
    my ($self, $SO_term) = @_;
    $self->{SO_term} = $SO_term if $SO_term;
    return $self->{SO_term};
}

sub predicate {
    my ($self, $predicate) = @_;
    $self->{predicate} = $predicate if $predicate;
    return $self->{predicate};
}

1;

package Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele;

use strict;
use warnings;

sub new {
    my ($class) = @;
    return bless {}, $class;
}

sub new_fast {
    my ($class, $hashref) = @;
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
       
        if (not $vfo->feature->isa('Bio::EnsEMBL::Transcript')) {
            $self->{affects_cds} = 0;
        }
        if (@{ $vfo->pep_coords } != 1) {
            $self->{affects_cds} = 0;
        }
        elsif @{ $vfo->pep_coords}[0]->isa('Bio::EnsEMBL::Mapper::Gap') {
            $self->{affects_cds} = 0;
        }
        else {

            my $tran = $vfo->feature;

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
