=head1 LICENSE

Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut


=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

 Questions may also be sent to the Ensembl help desk at
 <http://www.ensembl.org/Help/Contact>.

=cut
package Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::InitSubmitJob;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::BaseDataDumpsProcess');

sub fetch_input {
    my $self = shift;
	my @input;

	my $species   = $self->param('species');
	my $config    = $self->param('config');
    my $release   = $self->param('release');
    my $output_dir = $self->param('pipeline_dir');
    my $job_type  = $self->param('job_type'); # parse or dump

	my $script_args = {};
    my $input; 
    if ($job_type eq 'dump_population') {
        foreach my $dump_type (keys %$config) { # generic, sets, incl_consequences, svs, populations, individuals
            if ($dump_type eq 'populations') {
                my $prefetched_frequencies_dir = $self->param('prefetched_frequencies');
                my $input_gvf = "$output_dir/gvf/$species/$species.gvf.gz";
                foreach my $group (keys %{$config->{populations}}) {
                    my $frequencies_dir = "$prefetched_frequencies_dir/$group";
                    foreach my $population_name (keys %{$config->{populations}->{$group}}) {
                        my $file_name  = $config->{populations}->{$group}->{$population_name}->{file};
                        my $short_name = $config->{populations}->{$group}->{$population_name}->{short};
                        my $script_arg = "--input_gvf $input_gvf --population $population_name --frequencies_dir $frequencies_dir --short_name $short_name";
                        $script_args->{$script_arg} = $file_name;
                    }
                }
            }
        }
    } else {
	    my @arguments          = map {'--' . $_} @{$config->{generic}};
	    my $generic_script_arg = join(' ', @arguments);

        foreach my $dump_type (keys %$config) { # generic, sets, incl_consequences, svs, populations, individuals
            if ($dump_type eq 'failed') {
                next if ($job_type eq 'parse');
            }
            if ($dump_type eq 'populations') {
                next if ($job_type eq 'dump' || $job_type eq 'dump_population');
                foreach my $group (keys %{$config->{populations}}) {
                    foreach my $population_name (keys %{ $config->{populations}->{$group} }) {
                        my $file_name  = $config->{populations}->{$group}->{$population_name}->{file};
                        my $script_arg = "--population $population_name";
                        $script_args->{$script_arg} = $file_name;
                    }
                }
            }
            if ($dump_type eq 'sets') {
                foreach my $set_name (keys %{$config->{sets}}) {
                    my $script_arg = "--set_name $set_name $generic_script_arg";
                    my $file_name = "$species\_$set_name";
                    $script_args->{$script_arg} = $file_name;
                }
            } elsif ($dump_type eq 'individuals') {
                foreach my $individual_name (keys %{$config->{individuals}}) {
                    my $file_name = $species . '_' . $config->{individuals}->{$individual_name};
                    my $script_arg = "--individual $individual_name $generic_script_arg";
                    $script_args->{$script_arg} = $file_name;
                }
            } else {
                my @arguments = map {'--' . $_} @{$config->{$dump_type}};
                my $script_arg = join(' ', @arguments);
                my $file_name = "$species\_$dump_type";
                $script_args->{$script_arg} = $file_name;
            }
        }
    }
    if ($job_type eq 'dump' || $job_type eq 'dump_population') {  
        $input = $self->get_input_gvf_dumps($script_args); 
    } elsif ($job_type eq 'parse') {
        $input = $self->get_input_gvf2vcf($script_args); 
    } else {
        die "Job type must be parse or dump. $job_type is not recognised.";
    }

    $self->param('input_for_submit_job', $input); 
}


sub get_input_gvf_dumps {
    my $self = shift;
    my $script_args = shift;

	my $file_type       = 'gvf';
    my $script_dir      = $self->param('script_dir');
	my $script          = '/export/release/dump_gvf.pl';
    my $output_dir      = $self->param('pipeline_dir');
    my $connection_args = '--registry ' . $self->param('registry_file');
    my $species = $self->param('species');
    my @input = ();
    
    if ($self->param('job_type') eq 'dump_population') {
        foreach my $script_arg (keys %$script_args) {	
		    my $file_name = $script_args->{$script_arg};
            my $params = {};
            my $output_file = "--gvf_file $output_dir/$file_type/$species/$file_name.$file_type";
            my $err = "$output_dir/$file_type/$species/$file_name-$file_type.err";
            my $out = "$output_dir/$file_type/$species/$file_name-$file_type.out";
            $params->{'species'}          = $species;
            $params->{'script'}           = "$script_dir/$script";
            $params->{'connection_args'}  = $connection_args;
            $params->{'script_args'}      = $script_arg;
            $params->{'gvf_file'}         = $output_file;
            $params->{'err'}              = $err;
            $params->{'out'}              = $out;
        
            push @input, $params;
	    }
    } else {
        my $slice_load = $self->get_slice_load($species);    

        foreach my $script_arg (keys %$script_args) {	
		    my $file_name = $script_args->{$script_arg};
		    foreach my $slice_set (keys %$slice_load) {
			    my @seq_region_ids = @{$slice_load->{$slice_set}};
			    my $seq_region_file = '';
                my $seq_regions = '';
			    if (@seq_region_ids) {
                    # store seq_region ids in file
				    $seq_regions = $seq_region_ids[0] . '_' . $seq_region_ids[-1];
                    $seq_region_file = "$output_dir/$file_type/$species/$seq_regions.txt";
                    my $fh = FileHandle->new($seq_region_file, 'w');                                
                    foreach my $seq_region_id (@seq_region_ids) {
                        print $fh $seq_region_id, "\n";
                    }
                    $fh->close();
			    }
			    my $params = {};
			    my $seq_regions_arg = ($seq_regions ? "-$seq_regions" : '');
			    my $output_file = "--gvf_file $output_dir/$file_type/$species/$file_name" . "$seq_regions_arg.$file_type";
			    my $err = "$output_dir/$file_type/$species/$file_name" . "$seq_regions_arg.err";
			    my $out = "$output_dir/$file_type/$species/$file_name" . "$seq_regions_arg.out";
			    $params->{'species'}          = $species;
			    $params->{'script'}           = "$script_dir/$script";
			    $params->{'connection_args'}  = $connection_args;
			    $params->{'seq_region_file'}  = ($seq_region_file ? "--seq_region_file $seq_region_file" : '');
			    $params->{'script_args'}      = $script_arg;
			    $params->{'gvf_file'}         = $output_file;
			    $params->{'err'}              = $err;
			    $params->{'out'}              = $out;
			
			    push @input, $params;
		    }
	    }
    }
    return \@input;
}

sub get_input_gvf2vcf {
    my $self = shift;
    my $script_args = shift;

	my $file_type       = 'vcf';
    my $script_dir      = $self->param('script_dir');
	my $script          = '/misc/release/gvf2vcf.pl';
    my $output_dir      = $self->param('pipeline_dir');
    my $connection_args = '--registry ' . $self->param('registry_file');
    my $species = $self->param('species');
    my @input = ();

    # don't forget to parse Populations and Individuals

	foreach my $script_arg (keys %$script_args) {	
		my $file_name = $script_args->{$script_arg};
        my $params = {};
        my $err = "$output_dir/$file_type/$species/$file_name-$file_type.err";
        my $out = "$output_dir/$file_type/$species/$file_name-$file_type.out";
        $params->{'species'}          = $species;
        $params->{'script'}           = "$script_dir/$script";
        $params->{'connection_args'}  = $connection_args;
        $params->{'script_args'}      = $script_arg;
        $params->{'err'}              = $err;
        $params->{'out'}              = $out;
      
        my $gvf_file_name = "$file_name.gvf.gz";
        my $vcf_file_name = "$file_name.vcf";

        if ($file_name =~ m/generic/) {
            $gvf_file_name = "$species.gvf.gz";
            $vcf_file_name = "$species.vcf";
        }

        my $gvf_file = "$output_dir/gvf/$species/$gvf_file_name";
        die "GVF file $gvf_file not found" unless(-e $gvf_file);
        my $vcf_file = "$output_dir/vcf/$species/$vcf_file_name";	

        $params->{'gvf_file'} = "--gvf_file $gvf_file";
        $params->{'vcf_file'} = "--vcf_file $vcf_file";
        push @input, $params;
	}
    return \@input;
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
    $self->dataflow_output_id($dump_input_parameters, 1);
    return;
}

1;
