package Bio::EnsEMBL::Variation::Pipeline::DataDumps::Report;

use strict;
use JSON;
use FileHandle;

use base ('Bio::EnsEMBL::Variation::Pipeline::DataDumps::BaseDataDumpProcess');

sub fetch_input {
    my $self = shift;
	my @input;
	my $config_file = $self->param('config_file');
	my $working_dir = $self->param('data_dump_dir');
	my $file_type   = $self->param('file_type');
	my $report_name = $self->param('report_name');

	my $fh = FileHandle->new($config_file, 'r');
	my $config_text = <$fh>;
	my $config = decode_json($config_text);
	$fh->close();

	my $fh = FileHandle->new("$working_dir/$report_name", 'w');
	foreach my $species (keys %$config) {
		opendir(DIR, "$working_dir/$file_type/$species") or die $!;			
		while (my $file = readdir(DIR)) {
			next if ($file =~ m/^\./);
			print $fh $file, "\n";		
		}	
	}
	$fh->close();
	closedir(DIR);
}


1;
