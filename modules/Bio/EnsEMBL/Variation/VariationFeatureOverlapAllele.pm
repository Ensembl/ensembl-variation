package Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele;

use strict;
use warnings;

use Scalar::Util qw(weaken);

sub new {
    my ($class) = @_;
    return bless {}, $class;
}

sub new_fast {
    my ($class, $hashref) = @_;
    my $self = bless $hashref, $class;
    # avoid a memory leak, because the vfo also has a reference to us
    weaken $self->{variation_feature_overlap} if $self->{variation_feature_overlap};
    return $self;
}

sub variation_feature_overlap {
    my ($self, $variation_feature_overlap) = @_;
    if ($variation_feature_overlap) {
        $self->{variation_feature_overlap} = $variation_feature_overlap;
        # avoid a memory leak, because the vfo also has a reference to us
        weaken $self->{variation_feature_overlap}
    }
    return $self->{variation_feature_overlap};
}

sub seq {
    my ($self, $seq) = @_;
    $self->{seq} = $seq if $seq;
    return $self->{seq};
}


sub consequences {
    my ($self, @new_consequences) = @_;
    
    $self->{consequences} ||= [];

    if (@new_consequences) {
        my $consequences = $self->{consequences};
        push @$consequences, @new_consequences; 
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
            last if $consequence->is_definitive;
        }
    }
}

1;
