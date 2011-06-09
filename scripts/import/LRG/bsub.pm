#! perl -w

###
# A wrapper for the bsub
# 

use strict;
use warnings;

package bsub;

our @ISA = "Executable";

# Some default parameters
sub defaults {
    return {
        'executable' => 'bsub',
        'memory' => 4000
    }; 
}

# Permitted AUTOLOAD routines
sub permitted {
    my $self = shift;
    return [
        @{$self->SUPER::permitted()},
        'memory',
        'jobname',
        'logout',
        'logerr',
        'job'
    ];
}

# The parameters are built from the specified characteristics
sub parameters {
    my $self = shift;

    my $parameters = [
        '-K',
        '-J ' . $self->jobname(),
        '-o ' . $self->logout(),
        '-e ' . $self->logerr(),
        "-R \"select[mem>" . $self->memory() . "] rusage[mem=" . $self->memory() . "]\"",
        '-M ' . $self->memory(),
        $self->job->command()
    ];
    
    return $parameters;
    
}

# Execute submits the job to the farm. This requires that all the necessary parameters have been set 
sub execute {
    my $self = shift;
    
    die ("A job has not been set") unless (defined($self->job()));
    
    # This will halt execution until the job has finished
    return $self->SUPER::execute();
}

1;
