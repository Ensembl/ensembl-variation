package Bio::EnsEMBL::Variation::Pipeline::DataDumps::SummariseValidation;

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

	my $summary_file_name = '';
	if ($file_type eq 'gvf') {
		$summary_file_name = 'SummaryValidateGVF.txt';
	} elsif ($file_type eq 'vcf') {
		$summary_file_name = 'SummaryValidateVCF.txt';
	}

	my $fh = FileHandle->new($config_file, 'r');
	my $config_text = <$fh>;
	my $config = decode_json($config_text);
	$fh->close();
	my $failed_dumps = 0;

	my $fh = FileHandle->new("$working_dir/$summary_file_name", 'w');

	foreach my $species (keys %$config) {
		my $path = "$working_dir/$file_type/$species";
		opendir(DIR, $path) or die $!;			
		while (my $file = readdir(DIR)) {
			if ($file =~ m/^Validate(.)*out$/) {
				my $count = 0;
				if ($file_type eq 'gvf') {
					$count = `grep 'No errors in this file' $path/$file`;				
					if ($count == 0) {
						print $fh "GVF validator reports errors for $file\n";
						$failed_dumps++;
					}
				} elsif ($file_type eq 'vcf') {

				}
			}
		}	
		closedir(DIR);
	}
	$fh->close();
	$self->warning("For $failed_dumps dumps report errors.") if ($failed_dumps > 0);
}


1;
