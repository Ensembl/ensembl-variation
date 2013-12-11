package Bio::EnsEMBL::Variation::Pipeline::DataDumps::JoinDump;

use strict;
use warnings;

use FileHandle;

use base ('Bio::EnsEMBL::Variation::Pipeline::DataDumps::BaseDataDumpProcess');

sub run {
	my $self = shift;
	my @input = ();
	my $species   = $self->param('species');
	my $dump_dir  = $self->param('data_dump_dir');
	my $file_type = 'gvf';	
		
	my $working_dir = "$dump_dir/$file_type/$species/";
	my $files = $self->get_files($working_dir, $file_type);
	if ($self->contains_unjoined_files($files)) {		
		$self->join_gvf($working_dir, $files);
	}
	foreach my $type (keys %$files) {
		my $params = {};
		$params->{file_name}   = "$species-$type";
		$params->{working_dir} = $working_dir;
		$self->warning('JoinDump: ' . $params->{file_name} . ' ' . $params->{working_dir});	
		push @input, $params;
	}
	$self->param('input_for_validation', \@input);	
}

sub get_files {
	my $self = shift;
	my $working_dir = shift;
	my $file_type = shift;	
	opendir(DIR, $working_dir) or die $!; 
	my %files;
	my ($species, $range, $type);
	while (my $file = readdir(DIR)) {
    	next if ($file =~ m/^\./);
    	if ($file =~ m/\.$file_type/) {
			$file =~ s/\.$file_type//g;
			my @file_name_parts =  split('-', $file);
			if (scalar @file_name_parts == 3) {
    			($species, $range, $type) = @file_name_parts;
    			$files{$type}{$range} = 1;
			} else {
	    		($species, $type) = @file_name_parts;
    			$files{$type}{1} = 1;
			}
		}
	}
	closedir(DIR);
	return \%files;	
}

sub contains_unjoined_files {
	my $self = shift;
	my $files = shift;
	foreach my $type (keys %$files) {
		return ((scalar keys %{$files->{$type}}) > 1);
	}
	return 0;	
} 

sub join_gvf {
	my $self  = shift;
	my $working_dir = shift;
	my $files = shift;
	my $species = $self->param('species');
	my $tmp_dir = $self->param('tmp_dir');	

	# check there are the same number of sub_files for each dump_type
	my $first_type     = (keys %$files)[0];
	my $sub_file_count = scalar keys %{$files->{$first_type}};

	foreach my $type (keys %$files) {
    	my $file_count = scalar keys %{$files->{$type}};
		die "file count differs for type: $type" unless ($file_count == $sub_file_count);	
	}

	# first assemble header: 
	# from first file complete header, then only sequence_region information
	foreach my $type (keys %$files) {
    	my $fh = FileHandle->new("> $working_dir/$species-$type.gvf");
    	# header ---------------------------------------------------------------------
    	my @sub_files = keys %{$files->{$type}};
    	my $first_file_range = shift @sub_files;
    	my $tmp_header_fh = FileHandle->new("< $working_dir/$species-$first_file_range-$type.gvf");
    	while (<$tmp_header_fh>) {
        	chomp;
        	my $line = $_;
        	if ($line =~ m/^#/) {
            	print $fh $line, "\n";
        	} else {last;}
    	}
    	$tmp_header_fh->close();

    	foreach my $file_range (@sub_files) {
			$self->warning("<$working_dir/$species-$file_range-$type.gvf");		
        	$tmp_header_fh = FileHandle->new("< $working_dir/$species-$file_range-$type.gvf");
        	while (<$tmp_header_fh>) {
            	chomp;
            	my $line = $_;
            	last unless ($line =~ m/^#/);
            	if ($line =~ m/^##sequence-region/) {
                	print $fh $line, "\n";
            	}
        	}
        	$tmp_header_fh->close();
   	 	}
    	# header ---------------------------------------------------------------------

    	# body -----------------------------------------------------------------------
    	# refresh sub_files
    	@sub_files = keys %{$files->{$type}};
    	my $id_count = 1;
    	foreach my $file_range (@sub_files) {
        	$tmp_header_fh = FileHandle->new("< $working_dir/$species-$file_range-$type.gvf");
        	while (<$tmp_header_fh>) {
            	chomp;
            	my $line = $_;
            	next if ($line =~ m/^#/);
            	# do some corrections here:
            	my $gvf_line = get_gvf_line(\$line);
            	$gvf_line->{attributes}->{ID} = $id_count;

	            $line = join("\t", map {$gvf_line->{$_}} (
				'seq_id', 
				'source', 
				'type', 
				'start', 
				'end', 
				'score', 
				'strand', 
				'phase'));
            	my $attributes = join(";", map{"$_=$gvf_line->{attributes}->{$_}"} keys %{$gvf_line->{attributes}});
            	print $fh $line, "\t", $attributes, "\n";   
            	$id_count++;
        	}
        	$tmp_header_fh->close();
			system("gzip $working_dir/$species-$file_range-$type.gvf");
			system("mv $working_dir/$species-$file_range-$type.gvf.gz $tmp_dir");
    	}
    	# body -----------------------------------------------------------------------
    	$fh->close();
	}
}

sub get_gvf_line {
    my $line = shift; 
    my $gvf_line;
    my ($seq_id, $source, $type, $start, $end, $score, $strand, $phase, $attrib) = split(/\t/, $$line);
    $gvf_line->{seq_id} = $seq_id;
    $gvf_line->{source} = $source;
    $gvf_line->{type} = $type;
    $gvf_line->{start} = $start;
    $gvf_line->{end} = $end;
    $gvf_line->{score} = $score;
    $gvf_line->{strand} = $strand;
    $gvf_line->{phase} = $phase;

    my @attributes = split(';', $attrib);

    foreach my $attribute (@attributes) {
        my ($key, $value) = split('=', $attribute);
        if ($value) {
            $gvf_line->{attributes}->{$key} = $value;
        }
    }
    return $gvf_line;
}

sub write_output {
	my $self = shift;
	$self->dataflow_output_id($self->param('input_for_validation'), 1);
	return;
}



1;
