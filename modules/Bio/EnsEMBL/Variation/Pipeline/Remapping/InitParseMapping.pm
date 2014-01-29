package Bio::EnsEMBL::Variation::Pipeline::Remapping::InitParseMapping;

use strict;

use base ('Bio::EnsEMBL::Hive::Process');

sub write_output {
    my $self = shift;

	my $bam_files_dir   = $self->param('bam_files_dir');
	my $fasta_files_dir = $self->param('fasta_files_dir');

	my @input;
	my $file_number = 1;

	opendir(my $dh, $fasta_files_dir) or die "Not a directory $fasta_files_dir";
	my $max_file_count = scalar(grep {$_ =~ m/\.fa$/} readdir($dh));
	$self->warning($max_file_count);
	closedir($dh);

	while ($file_number <= $max_file_count) {
		my $params = {};
		die "Fasta file $fasta_files_dir/$file_number.fa" unless (-e "$fasta_files_dir/$file_number.fa");
		die "Bam file $bam_files_dir/$file_number.bam missing" unless (-e "$bam_files_dir/$file_number.bam");
		push @input, {
			'file_number' => $file_number,
			'fasta_file'  => "$fasta_files_dir/$file_number.fa",
			'bam_file'    => "$bam_files_dir/$file_number.bam",
		};
		$self->warning("$fasta_files_dir/$file_number.fa");
		$self->warning("$bam_files_dir/$file_number.bam");
		$file_number++;
	}
	$self->dataflow_output_id(\@input, 2);
	1;
}


1;
