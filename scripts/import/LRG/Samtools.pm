#! perl -w

###
# A wrapper for the Samtools program suite
# 

use strict;
use warnings;
use LRG::Executable;

package Samtools;

our @ISA = "Executable";

# Some default parameters
sub defaults {
    return {
        'executable' => 'samtools',
        'program' => '',
        'extra_parameters' => []
    }; 
}

sub permitted {
    my $self = shift;
    return [
        @{$self->SUPER::permitted()},
        'extra_parameters',
        'program'
    ];
}

sub parameters {
    my $self = shift;
 
    my $parameters = [$self->program(),@{$self->extra_parameters()}];
    return $parameters;
}

1;
