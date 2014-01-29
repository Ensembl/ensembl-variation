package Bio::EnsEMBL::Variation::Pipeline::Remapping::RunMapping;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Hive::Process');

use FileHandle;

sub fetch_input {
	my $self = shift;
	1;
}

sub run {
	my $self = shift;
	
	my $tool_dir      = $self->param('tool_dir'); 	
	my $bwa_dir       = $self->param('bwa_dir');
	my $samtools_dir  = $self->param('samtools_dir');
	my $fasta_file    = $self->param('fasta_file');
	my $sam_file      = $self->param('sam_file');
	my $file_number   = $self->param('file_number');
	my $bam_files_dir = $self->param('bam_files_dir');

	my $sam2bam_err = "$bam_files_dir/sam2bam.$file_number.err";
	my $sam2bam_out = "$bam_files_dir/sam2bam.$file_number.out";

	my $bam_file = $self->param('bam_file');
	my $err_file = $self->param('err_file');
	my $out_file = $self->param('out_file');

	my $new_assembly_fasta_file_dir = $self->param('new_assembly_fasta_file_dir');
	my $new_assembly_fasta_file_name = $self->param('new_assembly_fasta_file_name');
	my $look_up_file = "$new_assembly_fasta_file_dir/$new_assembly_fasta_file_name";
	
	#my $cmd = "$tool_dir/$bwa_dir/bwa mem -a $look_up_file $fasta_file";
	my $cmd = "$tool_dir/$bwa_dir/bwa mem -a $look_up_file $fasta_file";

	$self->run_cmd("$cmd 1>$sam_file 2>$err_file");

	$cmd = "gzip $sam_file";
	$self->run_cmd($cmd);
	
	$cmd = "$tool_dir/$samtools_dir/samtools view -uS $sam_file.gz | $tool_dir/$samtools_dir/samtools sort - $bam_files_dir/$file_number"; 
	$self->run_cmd("$cmd 1>$sam2bam_out 2>$sam2bam_err");

	return 1;
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
