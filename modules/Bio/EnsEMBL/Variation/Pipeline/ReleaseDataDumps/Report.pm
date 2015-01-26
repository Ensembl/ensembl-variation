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
package Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::Report;

use strict;
use JSON;
use FileHandle;

use base ('Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::BaseDataDumpsProcess');

sub fetch_input {
    my $self = shift;
	my $pipeline_dir = $self->param('pipeline_dir');
	my $file_type   = $self->param('file_type');
	my $report_name = $self->param('report_name');

    my $all_species = $self->get_all_species();

	my $fh = FileHandle->new("$pipeline_dir/$report_name", 'w');
	foreach my $species (keys %$all_species) {
		opendir(DIR, "$pipeline_dir/$file_type/$species") or die $!;			
		while (my $file = readdir(DIR)) {
			next if ($file =~ m/^\./);
			print $fh $file, "\n";		
		}	
	}
	$fh->close();
	closedir(DIR);
}


1;
