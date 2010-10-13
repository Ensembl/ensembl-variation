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

1;


