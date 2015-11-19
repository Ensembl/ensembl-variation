=head1 LICENSE
<<<<<<< HEAD

Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

=======
Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
>>>>>>> Speed up pipeline ENSVAR-283.
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

use constant GLOBAL_VF_COUNT => 1_000; # 5_000_000 1_000
use constant MAX_VF_LOAD => 100; # 2_000_000 100
use constant VF_PER_SLICE => 10_000; # 2_000_000 10_000
my $max_length = 5e7;
my $overlap = 0;

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
    my $global_vf_count = $self->get_global_vf_count();
    if ($global_vf_count > GLOBAL_VF_COUNT) {
      my $covered_seq_regions = $self->get_covered_seq_regions();
      my $vf_distributions = $self->get_vf_distributions($covered_seq_regions);
      $input = $self->get_input_gvf_dumps($script_args, $vf_distributions);
    } else {
      $input = $self->get_input_gvf_dumps($script_args);
    }
  } elsif ($job_type eq 'parse') {
    $input = $self->get_input_gvf2vcf($script_args); 
  } else {
    die "Job type must be parse or dump. $job_type is not recognised.";
  }

  $self->param('input_for_submit_job', $input); 
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

  my $gvf_dir = "$output_dir/gvf/$species/";
  opendir(DIR, $gvf_dir) or die $!;

  while (my $gvf_file = readdir(DIR)) {
    next if ($gvf_file =~ m/^\./);
    next if ($gvf_file =~ m/failed/); # don't parse gvf files storing failed variants
    if ($gvf_file =~ m/\.gvf\.gz$/) {
      my $script_arg = $self->get_script_arg($gvf_file, $script_args);
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
  closedir(DIR);
  return \@input;
}

sub get_script_arg {
  my $self = shift;
  my $file_name = shift;
  my $script_args = shift;
  while (my ($script_arg, $dump_type) = each %$script_args) {
    if ($file_name =~ /$dump_type/) {
      return $script_arg;
    }    
  }
}

sub get_input_gvf_dumps {
  my $self = shift;
  my $script_args = shift;
  my $vf_distributions = shift; 

  my $file_type       = 'gvf';
  my $script_dir      = $self->param('script_dir');
  my $script          = '/export/release/dump_gvf.pl';
  my $output_dir      = $self->param('pipeline_dir');
  my $connection_args = '--registry ' . $self->param('registry_file');
  my $species = $self->param('species');
  my @input = ();
  my $default_params = {
    'species' => $species,
    'script' => "$script_dir/$script",
    'connection_args' => $connection_args,
    'file_type' => $file_type,
    'output_dir' => $output_dir,
  };
  
  if ($vf_distributions) {
    foreach my $script_arg (keys %$script_args) {
      my $file_name = $script_args->{$script_arg};
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
  my $self = shift;
  my $covered_seq_regions_counts = shift;
  my $output_dir = $self->param('pipeline_dir');
  my $species = $self->param('species');
  my @vf_loads = ();

  my $current_vf_load = 0;
  my @seq_region_ids = ();

  while (my ($seq_region_id, $vf_count) = each %$covered_seq_regions_counts) {
    if ($vf_count > VF_PER_SLICE) {
      push @vf_loads, @{$self->get_split_slices($seq_region_id)};
    } else {
      if (($current_vf_load + $vf_count) > MAX_VF_LOAD) {
        push @seq_region_ids, $seq_region_id;
        push @vf_loads, $self->get_seq_regions(\@seq_region_ids, "$output_dir/gvf/$species/"); 
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
  }

  return \@vf_loads;
}

sub get_seq_regions {
  my $self = shift;
  my $seq_region_ids = shift;
  my $species_dir = shift;
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
  my $slice_pieces = split_Slices([$slice], $max_length, $overlap);
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

sub write_output { 
  my $self = shift;
  my $dump_input_parameters = $self->param('input_for_submit_job');
  $self->dataflow_output_id($dump_input_parameters, 1);
  return;
}

1;
