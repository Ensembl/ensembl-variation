#! perl -w

###
# A wrapper for the SMALT program
# 

use strict;
use warnings;
use LRG::Executable;

package SMALT;

our @ISA = "Executable";

# Some default parameters
sub defaults {
    return {
        'executable' => 'smalt',
        'extra_parameters' => ['-a']
    }; 
}

sub permitted {
    my $self = shift;
    return [
        @{$self->SUPER::permitted()},
        'hash',
        'inputfile',
        'outputfile',
        'matefile',
        'extra_parameters'
    ];
}

sub parameters {
    my $self = shift;
 
    my $parameters = ['map'];
    # Append any additional parameters specified
    push(@{$parameters},@{$self->extra_parameters()});
    
    # If an output file has been specified, add that as a parameter
    push(@{$parameters},'-o ' . $self->outputfile()) if (defined($self->outputfile()));
    
    # Lastly, append the hash and the input and, possibly, mate pair file
    push(@{$parameters},($self->hash(),$self->inputfile()));
    push(@{$parameters},$self->matefile()) if (defined($self->matefile()));
    
    return $parameters;
}

1;
