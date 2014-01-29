package Bio::EnsEMBL::Variation::Pipeline::Remapping::PreRunChecks;

use strict;
use warnings;

use FileHandle;
use Bio::EnsEMBL::Registry;

use base ('Bio::EnsEMBL::Hive::Process');


sub fetch_input {
	my $self = shift;
	1;	
}

sub run {
	my $self = shift;

	# 1. check all folders are created 
	my $pipeline_dir = $self->param('pipeline_dir');
	die "$pipeline_dir doesn't exist" unless (-d $pipeline_dir);		

	foreach my $folder qw/fasta_files_dir bam_files_dir old_assembly_fasta_file_dir new_assembly_fasta_file_dir mapping_results_dir/ {
		my $dir = $self->param($folder);
		die "$dir for $folder doesn't exist" unless (-d $dir);		
	}

	# 2. check new_assembly_fasta_file is indexed
	my $dir = $self->param('new_assembly_fasta_file_dir');
	foreach my $file_type (('.fa.amb', '.fa.ann', '.fa.bwt', '.fa.pac', '.fa.sa')) {
		unless ($self->count_files($dir, $file_type)) {
			die("New assembly file is not indexed. $file_type is missing.");
		}
	}
	# 3. bam_files_dir, mapping_results_dir empty?
	foreach my $name ('bam_files_dir') {
		$dir = $self->param($name);
		if (! $self->is_empty($dir)) {
			die("$name is not empty. Delete files before running the pipeline.");
		}
	}	
	# 4. if not generate_fasta_files check that fasta_files dir contains files
	$dir = $self->param('fasta_files_dir');
	unless ($self->param('generate_fasta_files')) {
		my $count = $self->count_files($dir, '.fa');
		if ($count == 0) {
			die ("There are no fasta_files. Set parameter 'generate_fasta_files' to 1 in the conf file.");
		}				
	}
	# 5. check bwa and samtools are working

	
	1;
}

sub is_empty {
	my $self = shift;
	my $dir = shift;
    opendir(my $dh, $dir) or die "Not a directory $dir";
    my $count =  scalar(grep { $_ ne "." && $_ ne ".." } readdir($dh)) == 0;
	closedir($dh);
	return $count;
}

sub count_files {
	my $self = shift;
	my $dir = shift;
	my $file_type = shift;
    opendir(my $dh, $dir) or die "Not a directory $dir";
    my $count = scalar(grep { $_ =~ m/\Q$file_type$/ } readdir($dh)) == 0;
	closedir($dh);
	return $count;
}


sub write_output {
	my $self = shift;
	1;
}



1;
