package Bio::EnsEMBL::Variation::Pipeline::DataDumps::InitSubmitJob;

use strict;
use warnings;
use JSON;
use FileHandle;
use File::Path qw(make_path);

use base ('Bio::EnsEMBL::Variation::Pipeline::DataDumps::BaseDataDumpProcess');

sub fetch_input {
    my $self = shift;
	my @input;

	my $species   = $self->param('species');
	my $config    = $self->param('config');
    my $release   = $self->param('release');

	$self->warning("Create gvf dump jobs for: $species");

	my $file_type = $self->param('file_type');
	my $job_type  = $self->param('job_type'); # dump or parse
	my $script = '';
	if ($job_type eq 'dump') {
		$script = '/export/dump_gvf_fast.pl';
	} elsif ($job_type eq 'parse') {
		$script = '/misc/gvf2vcf_08_20.pl';
	} else {
		die "Job type must be parse or dump. $job_type is not recognised.";
	}

    my $output_dir      = $self->param('data_dump_dir');
    my $connection_args = '--registry ' . $self->param('ensembl_registry');

	my $slice_load->{1} = []; # inititalise for parse jobs;
	if ($job_type eq 'dump') {
		$slice_load = $self->get_slice_load($species);    
	}

	my $script_args        = {};
	my @arguments          = map {'--' . $_} split(',', $config->{generic});
	my $generic_script_arg = join(' ', @arguments);

	foreach my $dump_type (keys %$config) {
		my $script_arg = '';
		if ($dump_type eq 'populations') {
			my $tmp_population_freq_dir = "$output_dir/prefetched_frequencies/1000G/";
			foreach my $desc (keys %{$config->{populations}}) {	
				my $file_name = $config->{populations}->{$desc};
				my @descs = split('_', $desc);
				my $short_name = $descs[2];
				$script_arg = "--population $desc $generic_script_arg --tmp_dir $tmp_population_freq_dir --cache_file $short_name\_frequencies_1000G_SEQNAME.txt ";
				$script_args->{$script_arg} = $file_name;
			}
		} elsif ($dump_type eq 'sets') {
			foreach my $desc (keys %{$config->{sets}}) {
				my $file_name = $config->{sets}->{$desc};
				$script_arg = "--set $desc $generic_script_arg";
				$script_args->{$script_arg} = $file_name;
			}
		} elsif ($dump_type eq 'individuals') {
			foreach my $desc (keys %{$config->{individuals}}) {
				my $file_name = $config->{individuals}->{$desc};
				$script_arg = "--individual $desc $generic_script_arg";
				$script_args->{$script_arg} = $file_name;
			}
		} else {
			my @arguments = map {'--' . $_} split(',', $config->{$dump_type});
			$script_arg = join(' ', @arguments);
			$script_args->{$script_arg} = $dump_type;
		}
	}

	foreach my $script_arg (keys %$script_args) {	
		my $file_name = $script_args->{$script_arg};
		foreach my $slice_set (keys %$slice_load) {
			my @seq_region_ids = @{$slice_load->{$slice_set}};
			my $seq_regions = '';
			if (@seq_region_ids) {
				# only provide the range: first and last, script figures out middle part
				$seq_regions = $seq_region_ids[0] . '_' . $seq_region_ids[-1];
			}
			my $params = {};
			my $seq_regions_arg = ($seq_regions ? "$seq_regions-" : '');
			my $output_file = "--output_file $output_dir/$file_type/$species/$species-" . $seq_regions_arg . "$file_name.$file_type";
			my $err = "$output_dir/$file_type/$species/$species-" . $seq_regions_arg . "$file_name.err";
			my $out = "$output_dir/$file_type/$species/$species-" . $seq_regions_arg . "$file_name.out";
			$params->{'seq_regions_arg'}  = $seq_regions_arg;
			$params->{'species'}          = $species;
			$params->{'script'}           = $self->param('script_dir') . $script;
			$params->{'connection_args'}  = $connection_args;
			$params->{'seq_region_range'} = ($seq_regions ? "--seq_region_range $seq_regions" : '');
			$params->{'script_args'}      = $script_arg;
			$params->{'output_file'}      = $output_file;
			$params->{'err'}              = $err;
			$params->{'out'}              = $out;
			
			if ($job_type eq 'parse') {
				my $gvf_file_name = '';
				my $vcf_file_name = '';
				my $uc_species = ucfirst $species;
				if ($file_name eq 'generic') {
					$gvf_file_name = "$uc_species.gvf.gz";
					$vcf_file_name = "$uc_species.vcf";
				} else {
					$gvf_file_name = "$uc_species\_$file_name.gvf.gz";
					$vcf_file_name = "$uc_species\_$file_name.vcf";
				}
				my $gvf_file = "$output_dir/gvf/$species/$gvf_file_name";
				die "GVF file $gvf_file not found" unless(-e $gvf_file);
				my $vcf_file = "$output_dir/vcf/$species/$vcf_file_name";	
				$params->{'gvf_file'} = "--gvf_file $gvf_file";
				$params->{'vcf_file'} = "--vcf_file $vcf_file";
				$params->{'output_file'} = '';
			}

			push @input, $params;
		}
	}
    $self->param('input_for_submit_job', \@input); 
}

sub get_slice_load {
	my $self = shift;
	my $species = shift;
    my $global_vfs_count = $self->_get_global_vf_count($species);
    my $count = 0;
    my $instance = 1;
	my $slice_load;
    if ($global_vfs_count > 15000000) {
        my $vf_counts = $self->_get_vf_counts_per_slice($species);
        my $max_load = $global_vfs_count / 20;
        my @seq_region_ids = @{$self->_get_seq_region_ids($species)};
        foreach my $seq_region_id (@seq_region_ids) {
            if (defined $vf_counts->{$seq_region_id}) {
                my $local_vfs_count = $vf_counts->{$seq_region_id};
                $count += $local_vfs_count;
                if ($count < $max_load) {
                    push @{$slice_load->{$instance}}, $seq_region_id;
                } else {
                    push @{$slice_load->{$instance}}, $seq_region_id;
                    $instance++;
                    $count = 0;
                }
            }
        }
	 } else {
        push @{$slice_load->{$instance}}, ();
    }
	return $slice_load;	
}

sub _get_seq_region_ids {
    my $self = shift;
    my $species = shift;
    my $cdba = $self->get_species_adaptor($species, 'core');
    my $sa = $cdba->get_SliceAdaptor;
    my $slices = $sa->fetch_all('toplevel', undef, 0, 1);
    my @seq_region_ids;
    foreach my $slice (@$slices) {
        push @seq_region_ids, $slice->get_seq_region_id;
    }
    my @sorted_ids = sort {$a <=> $b} @seq_region_ids;
    return \@sorted_ids;
}

sub _get_global_vf_count {
    my $self = shift;
    my $species = shift;
    my $count;
    my $vdba = $self->get_species_adaptor($species, 'variation');
    my $dbh = $vdba->dbc->db_handle;
    my $sth = $dbh->prepare(qq{ SELECT count(*) FROM variation_feature; });
    $sth->{'mysql_use_result'} = 1;
    $sth->execute();
    $sth->bind_columns(\$count);
    $sth->fetch(); 
    $sth->finish();
    return $count;
}

sub _get_vf_counts_per_slice {
    my $self = shift;
    my $species = shift;
    my $counts;
    my $vdba = $self->get_species_adaptor($species, 'variation');
    my $dbh = $vdba->dbc->db_handle;
    my $sth = $dbh->prepare(qq{
        SELECT sr.seq_region_id, count(*)
        FROM seq_region sr, variation_feature vf
        WHERE sr.seq_region_id = vf.seq_region_id
        GROUP BY sr.seq_region_id;
    });
    $sth->{'mysql_use_result'} = 1;
    $sth->execute();
    my ($slice_id, $count);
    $sth->bind_columns(\$slice_id, \$count);
    while ($sth->fetch()) {
        if ($count > 0) {
            $counts->{$slice_id} = $count;
        }
    }
    $sth->finish();
    return $counts;
}

# some dataflow to perform
sub write_output { 
    my $self = shift;
    my $dump_input_parameters = $self->param('input_for_submit_job');
    $self->warning(scalar(@$dump_input_parameters) . ' data dump jobs have been created');
    $self->dataflow_output_id($dump_input_parameters, 1);
    return;
}

1;
