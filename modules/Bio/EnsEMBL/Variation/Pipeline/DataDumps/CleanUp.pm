package Bio::EnsEMBL::Variation::Pipeline::DataDumps::CleanUp;

use strict;
use JSON;
use FileHandle;


use base ('Bio::EnsEMBL::Hive::Process');

sub fetch_input {
	my $self = shift;
}

sub run {
    my $self = shift;
	my $data_dump_dir = $self->param('data_dump_dir');
	my $tmp_dir = $self->param('tmp_dir');

	my $config_file = $self->param('config_file');
	my $fh = FileHandle->new($config_file, 'r');	
	my $config_text = <$fh>;
	my $config = decode_json($config_text);
	$fh->close();

	my $file_type = $self->param('file_type');
	
	foreach my $species (keys %$config) {
		my $working_dir = "$data_dump_dir/$file_type/$species/";			
		`mv -f $working_dir/*.out $tmp_dir`;
		`mv -f $working_dir/*.err $tmp_dir`;			
		opendir(DIR, $working_dir) or die $!;
		while (my $file = readdir(DIR))	{
			next if ($file =~ m/^\./);
			if ($file =~ m/gvf/) {
				my $current_name = $file;
				$file =~ s/\.gvf//;					
				my ($species, $type) = split('-', $file);
				$species = ucfirst $species;
				if ($type eq 'generic')	{
					$file = $species . '.gvf';
				} else {
					$file = $species . '_' . $type . '.gvf';
				}
	
				`mv $working_dir/$current_name $working_dir/$file`;
				`gzip $working_dir/$file`;
			} elsif ($file =~ m/vcf$/) {
				`mv -f $working_dir/$file $tmp_dir`;
			} 
		}				
		closedir(DIR);
	}	
};

sub write_output {
	my $self = shift;
}

1;


















