=head1 LICENSE

Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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
package Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::CleanUp;

use strict;

use base ('Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::BaseDataDumpsProcess');

sub fetch_input {
	my $self = shift;
}

sub run {
    my $self = shift;
	my $data_dump_dir = $self->param('pipeline_dir');
	my $tmp_dir = $self->param('tmp_dir');
	my $file_type = $self->param('file_type');
    my $all_species = $self->get_all_species();	

	foreach my $species (keys %$all_species) {
		my $working_dir = "$data_dump_dir/$file_type/$species/";			
		`mv -f $working_dir/*.out $tmp_dir`;
		`mv -f $working_dir/*.err $tmp_dir`;			
		opendir(DIR, $working_dir) or die $!;
		while (my $file = readdir(DIR))	{
			next if ($file =~ m/^\./);
			if ($file =~ m/gvf/) {
                if ($file =~ m/generic/) {
				    my $current_name = $file;
                    $file =~ s/_generic//;
				    `mv $working_dir/$current_name $working_dir/$file`; # rename file name
                }
				`gzip $working_dir/$file`;
			} elsif ($file =~ m/\.vcf$/) {
				`mv -f $working_dir/$file $tmp_dir`;
			} 
		}				
		closedir(DIR);
	}	
}

sub write_output {
	my $self = shift;
}

1;

