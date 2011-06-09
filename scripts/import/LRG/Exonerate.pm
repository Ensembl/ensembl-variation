#! perl -w

###
# A wrapper for the Exonerate program
# 

use strict;
use warnings;
use LRG::Executable;

package Exonerate;

our @ISA = "Executable";

# Some default parameters
sub defaults {
    return {
        'executable' => 'exonerate',
        'extra_parameters' => []
    }; 
}

sub permitted {
    my $self = shift;
    return [
        @{$self->SUPER::permitted()},
        'query',
        'target',
        'model',
        'extra_parameters'
    ];
}

sub parameters {
    my $self = shift;
 
    # Check that query and target have been specified
    die ("A query file must be given as input to Exonerate") unless (defined($self->query()));
    die ("A target file must be given as input to Exonerate") unless (defined($self->target()));
    
    my $parameters = ['--query',$self->query(),'--target',$self->target()];
    push(@{$parameters},('--model',$self->model())) if (defined($self->model()));
    push(@{$parameters},@{$self->extra_parameters()});
    
    return $parameters;
}

1;
