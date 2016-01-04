=head1 LICENSE

Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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
package Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::SpeciesFactory;

use strict;
use JSON;
use FileHandle;
use File::Path qw(make_path);

use base ('Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::BaseDataDumpsProcess');

sub run {
    my $self = shift;
	my @input;
	my $config_file = $self->param('config_file');
	my $fh = FileHandle->new($config_file, 'r');
	my $config_text = <$fh>;
	my $config = decode_json($config_text);
	$fh->close();	
	foreach my $species (keys %$config) {
		my $params = {};
		$params->{species} = $species;
		$params->{config}  = $config->{$species};
		push @input, $params;			
	}
    $self->param('species', \@input); 
}

# some dataflow to perform
sub write_output { 
    my $self = shift;
    my $input = $self->param('species');
    $self->dataflow_output_id($input, 2);
    return 1;
}

1;
