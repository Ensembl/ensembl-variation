package Bio::EnsEMBL::Variation::OverlapConsequence;

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
