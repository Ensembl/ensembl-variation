=head1 LICENSE
Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2021] EMBL-European Bioinformatics Institute
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

package Bio::EnsEMBL::Variation::Pipeline::Remapping::FilterMapping;

use strict;
use warnings;

use FileHandle;
use Bio::EnsEMBL::Registry;

use base ('Bio::EnsEMBL::Variation::Pipeline::Remapping::BaseRemapping');
use Bio::EnsEMBL::Variation::Utils::RemappingUtils qw(filter_read_mapping);

sub fetch_input {
  my $self = shift;
  # initialise file names
  my $params = {
    mapping_results_dir   => 'file_mappings',
    filtered_mappings_dir => 'file_filtered_mappings',
    statistics_dir        => 'file_statistics',
    dump_features_dir     => 'file_init_feature',
    load_features_dir     => 'file_load_features',
    fasta_files_dir       => 'fasta_file',
  };
  my $file_number = $self->param('file_number');
  my $individual_dir = '';
  if ($self->param('mode') eq 'remap_read_coverage') {
    my $individual_id = $self->param('individual_id');
    $individual_dir = "/$individual_id/";
  }
  foreach my $param (keys %$params) {
    my $dir = $self->param($param) . $individual_dir;
    if ($param =~ /mapping_results_dir/) {
      my $file_mappings = "$dir/mappings_$file_number.txt";
      my $file_failed_mappings = "$dir/failed_mapping_$file_number.txt";
      $self->param('file_mappings', $file_mappings);
      $self->param('file_failed_mappings', $file_failed_mappings);
    } elsif ($param =~ /filtered_mappings_dir/) {
      my $file_filtered_mappings = "$dir/$file_number.txt";
      my $file_failed_filtered_mappings = "$dir/failed_mappings_$file_number.txt";
      $self->param('file_filtered_mappings', $file_filtered_mappings);
      $self->param('file_failed_filtered_mappings', $file_failed_filtered_mappings);
    } elsif ($param =~ /fasta_files_dir/) {
      my $fasta_file = "$dir/$file_number.fa";
      $self->param('fasta_file', $fasta_file);
    } else {
      my $file = "$dir/$file_number.txt";
      my $file_param = $params->{$param};
      $self->param($file_param, $file);
    }
  }
}

sub run {
  my $self = shift;
}

sub write_output {
  my $self = shift;
}

sub get_seq_region_names {
  my $self = shift;
  my $coord = shift;
  my $cdba = $self->get_newasm_core_database_connection;
  my $seq_region_names = {};
  my $sa = $cdba->get_SliceAdaptor;
  foreach my $slice (@{$sa->fetch_all($coord)}) {
    $seq_region_names->{$slice->seq_region_name} = 1;
  } 
  return $seq_region_names; 
}

sub get_old_chrom_name {
  my $self = shift;
  my $query_name = shift;
  my @query_name_components_fst_part = split('-', $query_name, 5);
  my @query_name_components_snd_part = split(':', $query_name_components_fst_part[4], 2);
  return $query_name_components_snd_part[0];
}

sub write_statistics {
  my $self = shift;
  my $file_statistics = $self->param('file_statistics');
  my $fh_statistics = FileHandle->new($file_statistics, 'w');

  foreach my $stats (qw/count_input_ids pre_count_unmapped pre_count_mapped stats_failed stats_unique_map stats_multi_map stats_unique_map_after_filter stats_multi_map_after_filter stats_failed_poor_score stats_failed_after_filter stats_failed_non_chrom stats_exceeds_max_map_weight/) {
    if ($self->param($stats)) {
      my $count = $self->param($stats);
      print $fh_statistics "$stats=$count\n";
    }
  }
  $fh_statistics->close();
}

sub read_data {
  my $self = shift;
  my $line = shift;
  my @key_values = split("\t", $line);
  my $mapping = {};
  foreach my $key_value (@key_values) {
    my ($table_name, $value) = split('=', $key_value, 2);
    $mapping->{$table_name} = $value;
  }
  return $mapping;
}

sub print_feature_line {
  my ($self, $query_name, $new_seq_info, $score) = @_;
  my ($seq_name, $start, $end, $strand) = split(' ', $new_seq_info);
  my $line =  join("\t", ($query_name, $seq_name, $start, $end, $strand, $score));
  return $line;
}

# returns the row that will be loaded into the new variation_feature_table
# change variation_feature_id to old_variation_feature_id
sub print_complete_feature_line {
  my $self = shift;
  my $data = shift;
  my $variation_feature_id = $data->{variation_feature_id};
  $data->{variation_feature_id_old} = $variation_feature_id;
  my @output = ();
#  $self->warning(join(', ', keys %$data));
  foreach my $column_name (sort keys %$data) {
    unless ($column_name =~ /^variation_feature_id$/ || $column_name =~ /^seq_region_name$/) {
      if ($self->param('mode') eq 'remap_QTL') {
        unless ($column_name =~ /^entry$/ || $column_name =~ /^variation_feature_id_old$/) {
          push @output, $data->{$column_name};
        }
      } else {
        push @output, $data->{$column_name};
      }
    }
  }
  my $line = join("\t", @output);
  return $line;
}

sub best_mapping {
  my $self = shift;
  my $threshold = shift;
  my $prev_location = shift;
  my $new_mappings = shift;
  my $new_mapping = '';
  my $diff = 1_000_000_000;

  foreach my $start (keys %$new_mappings) {
    my $score = $new_mappings->{$start};
    if ($score >= $threshold) {
      my $new_diff = abs($prev_location - $start);
      if ($diff > $new_diff) {
        $diff = $new_diff;
        $new_mapping = $start;
      }
    }
  }
  return $new_mapping;
}

sub coords_are_in_order {
  my $self = shift;
  my $order = shift;
  my $mappings = shift;

  for (my $i = 0; $i < scalar @$order - 1; $i++) {
    my $a_coord = $order->[$i];
    my $b_coord = $order->[$i + 1];
    my $a_start = get_start($mappings->{$a_coord});
    my $b_start = get_start($mappings->{$b_coord});
    if ($a_start && $b_start) {
      if ($a_start > $b_start) {
        return 0;
      }
    }
  }
  return 1;
}

sub get_start {
  my $self = shift;
  my $coord = shift;
  if ($coord) {
    my @keys = keys %$coord;
    return $keys[0];
  }
  return 0;
}

sub overlap {
  my $self = shift;
  my $a = shift;
  my $b = shift;
  return (scalar @$a == scalar @$b);
}

sub report_failed_read_mappings {
  my $self = shift;
  my $count_mapped = 0;
  my $count_unmapped = 0;
  my $file_mappings = $self->param('file_mappings');
  my $file_failed_mappings = $self->param('file_failed_mappings');

  my $mapped = {};
  my $unmapped = {};
  my $fh_mappings = FileHandle->new($file_mappings, 'r');
  while (<$fh_mappings>) {
    chomp;
    my ($query_name, $seq_name, $start, $end, $strand, $map_weight, $score) = split("\t", $_);
    my ($id, $type) = split(':', $query_name);
    $mapped->{$id} = 1;
  }
  $fh_mappings->close();

  $count_mapped = scalar keys %$mapped;

  $self->param('pre_count_mapped', $count_mapped);
  $self->param('pre_count_unmapped', $count_unmapped);
  my $count_input_ids = 0;
  my $fasta_file = $self->param('fasta_file');
  my $fh_fasta_file = FileHandle->new($fasta_file, 'r');
  while (<$fh_fasta_file>) {
    chomp;
    if (/^>/) {
      $count_input_ids++;
    }
  }
  $fh_fasta_file->close();
  $self->param('count_input_ids', $count_input_ids);
}

sub filter_read_mapping_results {
  my $self = shift;
  my $feature_type_id = shift; # phenotype_feature_id
  $feature_type_id ||= 'entry';

  my $config = {
    feature_type_id => $feature_type_id,
    file_init_feature => $self->param('file_init_feature'),
    file_filtered_mappings =>  $self->param('file_filtered_mappings'),
    file_mappings => $self->param('file_mappings'),
  };
  my $stats = filter_read_mapping($config);

  $self->param('stats_unique_map', $stats->{unique_map});
  $self->param('stats_multi_map', $stats->{multi_map});
  $self->param('stats_failed', $stats->{failed});
}

sub get_new_seq_region_ids {
  my $self = shift;
  my $seq_region_ids = {};
  my $cdba = $self->get_newasm_core_database_connection;
  my $sa = $cdba->get_SliceAdaptor;
  my $slices = $sa->fetch_all('toplevel', undef, 1);
  foreach my $slice (@$slices) {
    $seq_region_ids->{$slice->seq_region_name} = $slice->get_seq_region_id;
  }
  return $seq_region_ids;
}

1;
