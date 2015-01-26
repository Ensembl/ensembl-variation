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
package Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::PreRunChecks;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::BaseDataDumpsProcess');

use File::Path qw(make_path);

sub fetch_input {
	my $self = shift;
}

sub run {
    my $self = shift;
	my $file_type  = $self->param('file_type');	

    if ($file_type eq 'gvf') {
        foreach my $dir (qw/pipeline_dir script_dir/) {
            die "$dir doesn't exist" unless (-d $self->param($dir));
        }
		my $gvf_validator = $self->param('gvf_validator');
		my $so_file = $self->param('so_file');
		die "gvf validator not defined" unless (defined $gvf_validator);
		die "wrong location for gvf validator" unless (-f $gvf_validator);
		die "so file not defined" unless (defined $so_file);
		die "wrong location for so file" unless (-f $so_file);
    } elsif ($file_type eq 'vcf') {
        my $location = `which vcf-validator`;
		die "vcf-validator command not found" if ($location =~ m/Command not found/);				
    } else {
        die "File type: $file_type is not recognised. It must be gvf or vcf.";
    }

    $self->create_species_dir_tree();
}

sub write_output {
	my $self = shift;
}

sub create_species_dir_tree {
    my $self = shift;
    my $pipeline_dir = $self->param('pipeline_dir');
    my $file_type = $self->param('file_type');
	my $dump_dir = "$pipeline_dir/$file_type/";
	make_path($dump_dir) unless (-d $dump_dir);	

    my $all_species = $self->get_all_species();

	foreach my $species (keys %$all_species) {
		if (-d "$dump_dir/$species") {
			unless (is_empty("$dump_dir/$species")) {
		        die("$dump_dir/$species is not empty. Delete files before running the pipeline.");
			}
		} else {
			make_path("$dump_dir/$species") or die "Failed to create dir $dump_dir/$species $!";
		}	
	}					
}

sub is_empty {
    my $dir = shift;
    opendir(my $dh, $dir) or die "Not a directory $dir";
    my $count =  scalar(grep { $_ ne "." && $_ ne ".." } readdir($dh)) == 0;
    closedir($dh);
    return $count;
}


1;
