=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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


=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

 Questions may also be sent to the Ensembl help desk at
 <http://www.ensembl.org/Help/Contact>.

=cut
package Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::SubmitJob;

use strict;
use FileHandle;
use base ('Bio::EnsEMBL::Hive::Process');

sub fetch_input {}

sub run {
    my $self = shift;
    my @args = ();
    my $script = $self->param('script');
    my $species = $self->param('species');
    push @args, "--species $species";
    my $debug = $self->param('debug') ? '--debug' : '';
    push @args, $debug;
    foreach my $arg (qw/connection_args seq_region_file script_args output_file gvf_file vcf_file/) {
        if (defined $self->param($arg)) {
            push @args, $self->param($arg);
        }
    }
    my $err = $self->param('err');
    my $out = $self->param('out');
	my $cmd = "perl $script " . join(' ', @args); 

	$self->run_cmd("$cmd 1>$out 2>$err");					
    return 1;
}

sub write_output {
    my $self = shift;
}

sub run_cmd {
	my $self = shift;
	my $cmd = shift;
	if (my $return_value = system($cmd)) {
		$return_value >>= 8;
		die "system($cmd) failed: $return_value";
	}
}


1;
