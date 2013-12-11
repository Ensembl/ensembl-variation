package Bio::EnsEMBL::Variation::Pipeline::DataDumps::ValidateDump;


use strict;
use base ('Bio::EnsEMBL::Hive::Process');

sub fetch_input {}

sub run {
	my $self = shift;

	my $gvf_validator = $self->param('gvf_validator');
	my $so_file       = $self->param('so_file');

	my $working_dir   = $self->param('working_dir');
	my $file_name     = $self->param('file_name');	
	
	my $err = "$working_dir/Validate\_$file_name.err";
	my $out = "$working_dir/Validate\_$file_name.out";

	my $file = "$working_dir/$file_name.gvf";
	my $cmd = "perl $gvf_validator --so_file $so_file $file";
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
