package Bio::EnsEMBL::Variation::Pipeline::DataDumps::PreRunChecks;

use strict;
use warnings;

use FileHandle;
use JSON;
use Bio::EnsEMBL::Registry;
use File::Path qw(make_path);

use base ('Bio::EnsEMBL::Hive::Process');


sub fetch_input {
	my $self = shift;
	1;	
}

=begin
Organise data dumps:
	-check directories exist
=end
=cut


sub run {
	my $self = shift;

	# directory: data_dumps/file_type/species/
	my $output_dir = $self->param('data_dump_dir');
	die "$output_dir doesn't exist" unless (-d $output_dir);
	my $tmp_dir    = $self->param('tmp_dir');
	die "$tmp_dir doesn't exist" unless (-d $tmp_dir);

	my $file_type  = $self->param('file_type');	

	if ($file_type eq 'gvf') {
		my $gvf_validator = $self->param('gvf_validator');
		my $so_file = $self->param('so_file');
		die "gvf validation tool not defined" unless ($gvf_validator);
		die "so file not defined" unless ($so_file);
	} elsif ($file_type eq 'vcf') {
		my $location = `which vcf-validator`;
		die "vcf-validator command not found" if ($location =~ m/Command not found/);				
	} else {
		die "File type: $file_type is not recognised. It must be gvf or vcf.";
	}
	
	# compare species in config and file_type dir
	my @species_from_config = @{$self->get_species_from_config()};

	my $dump_dir =  "$output_dir/$file_type/";
	make_path($dump_dir) unless (-d $dump_dir);	
	foreach my $species (@species_from_config) {
		if (-d "$dump_dir/$species") {
			if (!$self->is_empty("$dump_dir/$species")) {
				die("$dump_dir/$species is not empty. Delete files before running the pipeline.");
			}
		} else {
			make_path("$dump_dir/$species") or die "Failed to create dir $dump_dir/$species $!";
		}	
	}					
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

sub get_species_from_config {
	my $self = shift;
	my $config_file = $self->param('config_file');
	my $fh = FileHandle->new($config_file, 'r');
	my $config_text = <$fh>;
	$fh->close();		
	my $config = decode_json($config_text);
	my @species = (keys %$config);
	return \@species;
}

sub write_output {
	my $self = shift;
	1;
}


1;
