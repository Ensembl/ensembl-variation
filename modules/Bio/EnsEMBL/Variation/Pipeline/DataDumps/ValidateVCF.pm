package Bio::EnsEMBL::Variation::Pipeline::DataDumps::ValidateVCF;


use strict;
use base ('Bio::EnsEMBL::Hive::Process');
use File::Basename;

sub fetch_input {
	my $self = shift;
}

sub run {
	my $self = shift;

	my $vcf_file = $self->param('vcf_file');	
	$vcf_file =~ s/--vcf_file //;
	my ($file_name, $working_dir, $suffix) = fileparse($vcf_file, qr/\.[^.]*/);	

	my $err = "$working_dir/Validate\_vcf\_$file_name.err";
	my $out = "$working_dir/Validate\_vcf\_$file_name.out";

	my $cmd;
	# sort and bgzip
	$cmd = "vcf-sort < $vcf_file | bgzip > $vcf_file.gz";
	$self->run_cmd($cmd);
	# validate 
	$cmd = "vcf-validator $vcf_file.gz";
	$self->run_cmd("$cmd 1>$out 2>$err");	
}

sub run_cmd {
	my $self = shift;
	my $cmd = shift;
	if (my $return_value = system($cmd)) {
		$return_value >>=8;
		die "system($cmd) failed: $return_value";
	}
}


1;
