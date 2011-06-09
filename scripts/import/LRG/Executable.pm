#! perl -w

###
# A framework for wrapping around an executable, can also submit to the farm
# 

use strict;
use warnings;

package Executable;

#ÊUse autoload for get/set methods
our $AUTOLOAD;

sub DESTROY { }

# Some default values
sub defaults {
    return {
        'tmpdir' => '/tmp/'
    };
}

# Permitted AUTOLOAD routines
sub permitted {
    return [
        'executable',
        'tmpdir',
        'parameters',
        'pid'
    ];
}

sub AUTOLOAD {
    my $self = shift;
    my $type = ref($self) or die("$self is not an object");

    #ÊThe name of the called subroutine
    my $name = $AUTOLOAD;
    # Strip away the pre-pended package info
    $name =~ s/.*://;

    # Check that the subroutine should exist for this module
    unless (grep {/^$name$/} @{$self->_get_set('_permitted')}) {
        die "Can't access `$name' field in class $type";
    }

    # Call the get/set method
    return $self->_get_set('_' . $name,@_);
}

#ÊConstructor
sub new {
    
    my $class = shift;
    
    my $self = bless({},$class);
    $self->initialize(@_);
    
    return $self;
    
}

# Initialize some variables
sub initialize {
    my $self = shift;
    my %vals = @_;
    
    #ÊThe get/set methods that exist in this module
    my $permitted = $self->permitted();
    $self->_get_set('_permitted',$permitted);

    # Set the pid of this process
    $self->pid($$);

    #ÊSet the default values
    my %DEFAULTS = %{$self->defaults()};
    map { $self->$_($DEFAULTS{$_}) } keys(%DEFAULTS);
    
    #ÊSet any field passed in via the parameter hash
    map { $self->$_($vals{$_}) } keys(%vals);
     
}

# Internal get/set method
sub _get_set {
    my $self = shift;
    my $attribute = shift;
    my $value = shift;
    
    if (defined($value)) {
        $self->{$attribute} = $value;
    }
    
    return $self->{$attribute};
}

# Construct the command
sub command {
    my $self = shift;
    
    my $cmd = $self->executable();
    
    # Add the parameters
    map {$cmd .= " $_"} @{$self->parameters()};
    
    return $cmd;
    
}

# Execute
sub execute {
    my $self = shift;
    
    my $cmd = $self->command();
    
    # Should handle this better with capturing output etc.
    my $output = `$cmd`;
    
    return $output;
}

1;