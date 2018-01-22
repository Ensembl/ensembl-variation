=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

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
