=head1 LICENSE
Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute
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

use FileHandle;
use Bio::EnsEMBL::Utils::Slice qw(split_Slices);
use base ('Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::BaseDataDumpsProcess');

my $global_vf_count_in_species; #= 5_000_000; # if number of vf in a species exceeds this we need to split up dumps
my $max_vf_load; #= 2_000_000; # group slices together until the vf count exceeds max_vf_load
my $vf_per_slice; #= 2_000_000; # if number of vf exceeds this we split the slice and dump for each split slice
my $max_split_slice_length; #= 500_000;

my $overlap = 1;
my $debug = 0;

my $debug_fh;

sub fetch_input {
  my $self = shift;
  my @input;

  my $species   = $self->param('species');
  my $config    = $self->param('config');
  my $release   = $self->param('release');
  my $job_type  = $self->param('job_type'); # parse or dump
  $debug = $self->param('debug');
  my $output_dir = $self->data_dir($species);

  $global_vf_count_in_species = $self->param('global_vf_count_in_species') || $global_vf_count_in_species;
  $max_vf_load = $self->param('max_vf_load') || $max_vf_load; # group slices together until the vf count exceeds max_vf_load
  $vf_per_slice = $self->param('vf_per_slice') || $vf_per_slice; # if number of vf exceeds this we split the slice and dump for each split slice
  $max_split_slice_length = $self->param('max_split_slice_length') || $max_split_slice_length;

  my $fh;
  if ($debug) {
    $debug_fh =  FileHandle->new("$output_dir/$species\_initSubmitJob.txt", 'w');
  }

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
    my $global_vf_count = $self->get_global_vf_count();
    if ($debug) {
      print $debug_fh "GLOBAL_VF_COUNT\t$global_vf_count\n";
    }
    if ($global_vf_count > $global_vf_count_in_species) {
      my $covered_seq_regions = $self->get_covered_seq_regions();
#      if ($debug) {
#        foreach my $key (keys %$covered_seq_regions) {
#          print $debug_fh "COVERED_SE"$key, ' ', $covered_seq_regions->{$key}, "\n";
#        }   
#      }
      my $vf_distributions = $self->get_vf_distributions($covered_seq_regions,$species,$output_dir);
#      if ($debug) {
#        foreach my $distribution (@$vf_distributions) {
#          foreach my $key (keys %$distribution) {
#            print $fh $key, ' ', $distribution->{$key}, "\n";
#          }
#          print $fh "\n";
#        }
#      }
      $input = $self->get_input_gvf_dumps($script_args,$species,$output_dir,$vf_distributions);
    } else {
      $input = $self->get_input_gvf_dumps($script_args,$species,$output_dir);
    }
  } elsif ($job_type eq 'parse') {
    $input = $self->get_input_gvf2vcf($script_args,$species,$output_dir);
  } else {
    die "Job type must be parse or dump. $job_type is not recognised.";
  }
#  if ($debug) {
#    $debug_fh->close;
#    die "Finished writing report";
#  }
  $self->param('input_for_submit_job', $input); 
}

sub get_input_gvf2vcf {
  my ($self,$script_args,$species,$output_dir) = @_;

  my $file_type       = 'vcf';
  my $script_dir      = $self->param('script_dir');
  my $script          = '/misc/release/gvf2vcf.pl';
  my $connection_args = '--registry ' . $self->param('registry');
  my @input = ();

  # don't forget to parse Populations and Individuals

  my $gvf_dir = "$output_dir/gvf/$species/";
  opendir(my $dh, $gvf_dir) or die $!;
  my @dir_content = readdir($dh);
  closedir($dh);
  foreach my $gvf_file (@dir_content) {
    next if ($gvf_file =~ m/^\./);
    next if ($gvf_file =~ m/failed/); # don't parse gvf files storing failed variants
    if ($gvf_file =~ m/\.gvf\.gz$/) {
      my $script_arg = $self->get_script_arg($gvf_file, $script_args);
      $self->warning("$gvf_file $script_arg");
      my $params = {};
      my $file_name = $gvf_file;
      $file_name =~ s/\.gvf\.gz//;

      my $err = "$output_dir/$file_type/$species/$file_name.err";
      my $out = "$output_dir/$file_type/$species/$file_name.out";
      my $vcf_file = "$output_dir/$file_type/$species/$file_name.vcf";
      $params->{'species'}          = $species;
      $params->{'script'}           = "$script_dir/$script";
      $params->{'connection_args'}  = $connection_args;
      $params->{'script_args'}      = $script_arg;
      $params->{'err'}              = $err;
      $params->{'out'}              = $out;
      $params->{'gvf_file'}         = "--gvf_file $gvf_dir/$gvf_file";
      $params->{'vcf_file'}         = "--vcf_file $vcf_file";
      push @input, $params;
    }
  }
  return \@input;
}

sub get_script_arg {
  my ($self, $file_name, $script_args) = @_;
  my $return_script_arg = '';
  while (my ($script_arg, $dump_type) = each %$script_args) {
    $self->warning("get_script_arg $script_arg $dump_type");
    if ($file_name =~ /$dump_type/) {
      $return_script_arg =  $script_arg;
    }    
  }
  if ($return_script_arg) {
    return $return_script_arg;
  } else {
    die "Could not find script arg for $file_name"; 
  }

}

sub get_input_gvf_dumps {
  my ($self,$script_args,$species,$output_dir,$vf_distributions) = @_;

  my $file_type       = 'gvf';
  my $script_dir      = $self->param('script_dir');
  my $script          = '/export/release/dump_gvf.pl';
  my $connection_args = '--registry ' . $self->param('registry');

  my @input = ();
  my $run_in_debug_mode = $debug ? '--debug' : '';
  my $default_params = {
    'species' => $species,
    'script' => "$script_dir/$script",
    'connection_args' => $connection_args,
    'file_type' => $file_type,
    'output_dir' => $output_dir,
    'debug' => $run_in_debug_mode,
  };
  
  if ($vf_distributions) {
    foreach my $script_arg (keys %$script_args) {
      my $file_name = $script_args->{$script_arg};
      $self->warning("get_input_gvf_dumps $file_name $script_arg");
      foreach my $vf_distribution (@$vf_distributions) {   
        my $params = {};
        my $file_id = $vf_distribution->{file_id};
        my $output_file = "--$file_type\_file $output_dir/$file_type/$species/$file_name-$file_id.$file_type";
        my $err = "$output_dir/$file_type/$species/$file_name-$file_id.err";
        my $out = "$output_dir/$file_type/$species/$file_name-$file_id.out";
        $params->{'script_args'} = $script_arg;
        foreach my $param (qw/species script connection_args/) {
          $params->{$param} = $default_params->{$param};
        }
        if ($vf_distribution->{is_slice_piece}) {
          foreach my $param (qw/seq_region_id slice_piece_name slice_piece_start slice_piece_end/) {
            $params->{$param} = $vf_distribution->{$param};
          }
          $params->{is_slice_piece} = '--is_slice_piece';
        } else {
          $params->{seq_region_ids_file} = $vf_distribution->{seq_region_ids_file};
        }
        $params->{'gvf_file'} = $output_file;
        $params->{'err'} = $err;
        $params->{'out'} = $out;
        push @input, $params;
      }
    }
  } else {
    foreach my $script_arg (keys %$script_args) {
      my $params = {};
      my $file_name = $script_args->{$script_arg};
      my $file_id = $vf_distributions->{file_id};
      my $output_file = "--$file_type\_file $output_dir/$file_type/$species/$file_name.$file_type";
      $params->{'script_args'} = $script_arg;
      my $err = "$output_dir/$file_type/$species/$file_name.err";
      my $out = "$output_dir/$file_type/$species/$file_name.out";
      foreach my $param (qw/species script connection_args/) {
        $params->{$param} = $default_params->{$param};
      }
      $params->{'gvf_file'} = $output_file;
      $params->{'err'} = $err;
      $params->{'out'} = $out;
      push @input, $params;
    }
  }

  return \@input;
}

sub get_vf_distributions {
  my ($self, $covered_seq_regions_counts,$species,$output_dir) = @_;
  my @vf_loads = ();

  my $current_vf_load = 0;
  my @seq_region_ids = ();

  while (my ($seq_region_id, $vf_count) = each %$covered_seq_regions_counts) {
    if ($debug) {
      print $debug_fh "VF_PER_SLICE\t$seq_region_id\t$vf_count\n";
    }
    if ($vf_count > $vf_per_slice) {
      my @split_slices = @{$self->get_split_slices($seq_region_id)};
      if ($debug) {
        print $debug_fh "SPLIT_SLICES\t$seq_region_id\t", scalar @split_slices, "\n";
      }
      push @vf_loads, @split_slices;
    } else {
      if (($current_vf_load + $vf_count) > $max_vf_load) {
        push @seq_region_ids, $seq_region_id;
        push @vf_loads, $self->get_seq_regions(\@seq_region_ids, "$output_dir/gvf/$species/"); 
        if ($debug) {
          my $tmp_load = $current_vf_load + $vf_count; 
          print $debug_fh "JOIN_SLICES\t", join(',', @seq_region_ids), "\t$tmp_load\n";  
        }
        @seq_region_ids = ();
        $current_vf_load = 0;
      } else {
        push @seq_region_ids, $seq_region_id;
        $current_vf_load += $vf_count;
      }
    }
  }

  if (scalar @seq_region_ids > 0) {
    push @vf_loads, $self->get_seq_regions(\@seq_region_ids, "$output_dir/gvf/$species/"); 
    if ($debug) {
      print $debug_fh "JOIN_SLICES\t", join(',', @seq_region_ids), "\t$current_vf_load\n";  
    }

  }

  return \@vf_loads;
}

sub get_seq_regions {
  my ($self, $seq_region_ids, $species_dir) = @_;
  my $seq_regions = $seq_region_ids->[0] . '_' . $seq_region_ids->[-1];
  my $seq_region_ids_file = "$species_dir/$seq_regions.txt";
  my $fh = FileHandle->new($seq_region_ids_file, 'w');
  foreach my $seq_region_id (@$seq_region_ids) {
    print $fh $seq_region_id, "\n";
  }
  $fh->close();

  return {
    file_id => $seq_regions,
    seq_region_ids_file => "--seq_region_ids_file $seq_region_ids_file",
    is_slice_split => 0,
  };
}

sub get_split_slices {
  my $self = shift;
  my $seq_region_id = shift;
  my $species = $self->param('species');
  my $cdba = $self->get_species_adaptor($species, 'core');
  my $sa = $cdba->get_SliceAdaptor;
  my @slice_piece_ids = ();
  my $slice = $sa->fetch_by_seq_region_id($seq_region_id);
  my $seq_region_name = $slice->seq_region_name;
  my $slice_pieces = split_Slices([$slice], $max_split_slice_length, $overlap);
  foreach my $slice_piece (@$slice_pieces) {
    my $start = $slice_piece->start;
    my $end = $slice_piece->end;
    push @slice_piece_ids, {
      'is_slice_piece' => 1,
      'file_id' => "$seq_region_id\_$start\_$end",
      'seq_region_id' => "--seq_region_id $seq_region_id",
      'slice_piece_name' => "--seq_region_name $seq_region_name",
      'slice_piece_start' => "--slice_piece_start $start",
      'slice_piece_end' => "--slice_piece_end $end",
    };
  }
  return \@slice_piece_ids;
}

sub get_global_vf_count {
  my $self = shift;
  my $species = $self->param('species');
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

sub get_covered_seq_regions {
  my $self = shift;
  my $species = $self->param('species');
  my $counts;
  my $vdba = $self->get_species_adaptor($species, 'variation');
  my $cdba = $self->get_species_adaptor($species, 'core');
  my $toplevel_seq_region_ids = {};
  my $sa = $cdba->get_SliceAdaptor;
  my $toplevel_slices = $sa->fetch_all('toplevel');
  foreach my $toplevel_slice (@$toplevel_slices) {
    $toplevel_seq_region_ids->{$toplevel_slice->get_seq_region_id} = 1;
  }

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
      if ($toplevel_seq_region_ids->{$slice_id}) {
        $counts->{$slice_id} = $count;
      }
    }
  }
  $sth->finish();

  return $counts;
}

sub write_output { 
  my $self = shift;
  my $dump_input_parameters = $self->param('input_for_submit_job');
  $self->dataflow_output_id($dump_input_parameters, 2);
  return;
}

1;
