use strict;
use warnings;

package GVFDumper;

use base 'BaseDumper';

sub new {
	my $caller = shift;
	my $class = ref($caller) || $caller;
	my $self = $class->SUPER::new(@_);
	return bless $self, $class;
}

sub dump {
	my $self = shift;
	$self->SUPER::init();
    my %header = (no_sequence_region => 1);
    $self->SUPER::print_gvf_header(\%header);
    my $sub_slice;
    my $vf;
    my $svf;
	while ($sub_slice = $self->SUPER::next_slice()) {
        #print "DEBUG ", $sub_slice->name, " ", $sub_slice->start, " ", $sub_slice->end, "\n";
        unless ($self->{'just_svs'}) {
            my $vfs = $self->SUPER::vfs_on_slice($sub_slice);
            VF : while ($vf = $vfs->next()) {
                $self->SUPER::print_variation($vf);
            }
        }
        if ($self->{'include_svs'}) {
            my $svfs = $self->SUPER::svfs_on_slice($sub_slice);
            while ($svf = $svfs->next()) {
                 $self->SUPER::print_structural_variation($svf); 
            }
        }
	}
}
1;
