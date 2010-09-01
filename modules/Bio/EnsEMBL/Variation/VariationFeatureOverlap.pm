package Bio::EnsEMBL::Variation::VariationFeatureOverlap;

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

sub alt_alleles {
    my ($self, @new_alt_alleles) = @_;

    if (@new_alt_alleles) {
        my $alt_alleles = $self->{alt_alleles} ||= [];
        push @$alt_alleles, @new_alt_alleles;
    }

    return $self->{alt_alleles};
}

sub alleles {
    my ($self) = @_;

    my $alleles = $self->alt_alleles || [];
    
    unshift @$alleles, $self->reference_allele;

    return $alleles;
}

1;


