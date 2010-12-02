package Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
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
        weaken $self->{variation_feature_overlap};
    }
    
    return $self->{variation_feature_overlap};
}

sub feature_seq {
    # the sequence of this allele relative to the feature
    my ($self, $feature_seq) = @_;
    
    $self->{feature_seq} = $feature_seq if $feature_seq;
    
    unless ($self->{feature_seq}) {
        
        # check if we need to reverse complement the vf_seq
        
        my $vf = $self->variation_feature_overlap->variation_feature;
        my $feature = $self->variation_feature_overlap->feature;
        
        if ($vf->strand != $feature->strand) {
            my $vf_seq = $self->{variation_feature_seq};
            reverse_comp(\$vf_seq);
            $self->{feature_seq} = $vf_seq;
        }
        else {
            $self->{feature_seq} = $self->{variation_feature_seq};
        }
    }
    
    return $self->{feature_seq};
}

sub variation_feature_seq {
    # the sequence of this allele relative to the variation feature
    my ($self, $variation_feature_seq) = @_;
    $self->{variation_feature_seq} = $variation_feature_seq if $variation_feature_seq;
    return $self->{variation_feature_seq};
}

sub is_reference {
    my ($self, $is_reference) = @_;
    $self->{is_reference} = $is_reference if defined $is_reference;
    return $self->{is_reference};
}

sub dbID {
    my ($self, $dbID) = @_;
    $self->{dbID} = $dbID if defined $dbID;
    return $self->{dbID};
}

1;
