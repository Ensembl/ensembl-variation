package Bio::EnsEMBL::Variation::Pipeline::DataDumps::SpeciesFactory;

use strict;
use JSON;
use FileHandle;
use File::Path qw(make_path);

use base ('Bio::EnsEMBL::Variation::Pipeline::DataDumps::BaseDataDumpProcess');

sub run {
    my $self = shift;
	my @input;
	my $config_file = $self->param('config_file');
	my $fh = FileHandle->new($config_file, 'r');
	my $config_text = <$fh>;
	my $config = decode_json($config_text);
	$fh->close();	
	foreach my $species (keys %$config) {
		my $params = {};
		$params->{species} = $species;
		$params->{config} = $config->{$species};
		$self->warning('SpeciesFactory: ' . $species);	
		push @input, $params;			
	}
    $self->param('species', \@input); 
}

# some dataflow to perform
sub write_output { 
    my $self = shift;
    my $input = $self->param('species');
    $self->warning('Generate dumps for ' . scalar(@$input) . ' species');
    $self->dataflow_output_id($input, 2);
    return 1;
}

1;
