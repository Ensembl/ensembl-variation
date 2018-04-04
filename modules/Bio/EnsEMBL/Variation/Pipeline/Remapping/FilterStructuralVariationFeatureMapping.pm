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

package Bio::EnsEMBL::Variation::Pipeline::Remapping::FilterStructuralVariationFeatureMapping;

use strict;
use warnings;

use FileHandle;
use Bio::DB::Fasta;
use Bio::EnsEMBL::Registry;

use base ('Bio::EnsEMBL::Variation::Pipeline::Remapping::FilterMapping');

sub fetch_input {
  my $self = shift;
  $self->SUPER::fetch_input;
}

sub run {
  my $self = shift;
  $self->filter_svf_mapping_results();
  $self->join_svf_data();
  $self->SUPER::write_statistics();
}

sub write_output {
  my $self = shift;
}

sub filter_svf_mapping_results {
  my $self = shift;
  my $algn_score_threshold = $self->param('algn_score_threshold');

  my $file_mappings = $self->param('file_mappings');
  my $fh_mappings = FileHandle->new($file_mappings, 'r');

  my $lookup = {};
  my $file_number = $self->param('file_number');
  my $dump_features_dir = $self->param('dump_features_dir'); 
  my $fh_lookup = FileHandle->new("$dump_features_dir/lookup_$file_number.txt", 'r');
  while (<$fh_lookup>) {
    chomp;
    #93742   seq_region_name=X;seq_region_start=44601505;inner_start=44601505;inner_end=44619579;seq_region_end=44619579 
    my ($vf_id, $coords) = split/\t/;
    foreach my $coord (split(';', $coords)) {
      my ($key, $value) = split('=', $coord);
      $lookup->{$vf_id}->{$key} = $value;
    }
  }
  $fh_lookup->close();

  my ($stats_failed, $stats_unique_map, $stats_multi_map);

  my $mapped = {};
  while (<$fh_mappings>) {
    chomp;

    #feature_id-coord, coord: outer_start seq_region_start inner_start inner_end seq_region_end outer_end

    #107100:upstream 1       460501  461001  1       6       1       496     0.99001996007984   
    #query_name seq_region_name start end strand map_weight edit_dist score_count algn_score
    #102027-seq_region_start 10      44974586        44974758        1       145     1       274     score=0.850746268656716 173M28H
    my ($query_name, $new_seq_name, $new_start, $new_end, $new_strand, $map_weight, $algn_score, $cigar_str) = split("\t", $_); 

    my ($id, $coord) = split('-', $query_name);   
    if ($lookup->{$id}->{seq_region_name} eq "$new_seq_name") {
      $mapped->{$id}->{$coord}->{$new_start} = $algn_score;
    }
  }
  $fh_mappings->close();

  my $filtered_mappings = {};
  foreach my $id (keys %$mapped) {
    foreach my $coord_type (keys %{$mapped->{$id}}) {
      my $mappings = $mapped->{$id}->{$coord_type};
      my @coord_2_scores = sort { $mappings->{$b} <=> $mappings->{$a} } keys %$mappings;
      my $algn_score_threshold = 0.75;
      if ($mappings->{$coord_2_scores[0]} == 1.0) {
        $algn_score_threshold = 1.0;
      }
      my $count_exceed_threshold = grep {$_ >= $algn_score_threshold} values %$mappings;

      if ($count_exceed_threshold > 1) {
        my $prev_mapping = $lookup->{$id}->{$coord_type};
        my $mapped_location = $self->best_mapping($algn_score_threshold, $prev_mapping, $mappings);
        $filtered_mappings->{$id}->{$coord_type}->{$mapped_location} = $mappings->{$mapped_location};
      } else {
        foreach my $coord (@coord_2_scores) {
          my $score = $mappings->{$coord};
          if ($score >= $algn_score_threshold) {
            $filtered_mappings->{$id}->{$coord_type}->{$coord} = $score;
          }
        }
      }
    }
  }

  my $final_mappings = {};
  my $failed_mappings = {};
  foreach my $id (keys %$lookup) {
    my @prev_coord_types = grep {$_ ne 'seq_region_name'} keys %{$lookup->{$id}};
    my @new_coord_types =  keys %{$filtered_mappings->{$id}};
    if ($self->overlap(\@prev_coord_types, \@new_coord_types)) {
      my @start_order = qw/outer_start seq_region_start inner_start/;
      my @end_order = qw/inner_end seq_region_end outer_end/;
      my $start_coords_in_order = $self->coords_are_in_order(\@start_order, $filtered_mappings->{$id});
      my $end_coords_in_order = $self->coords_are_in_order(\@end_order, $filtered_mappings->{$id});
      if (!($start_coords_in_order && $end_coords_in_order)) {
        $failed_mappings->{$id} = 'Coords not in order';
        next;
      }
      # check start coords are smaller that end coords
      my $start = $self->get_start($filtered_mappings->{$id}->{seq_region_start});
      my $end = $self->get_start($filtered_mappings->{$id}->{seq_region_end});
      if ($end < $start) {
        my $swap_map = {
          'outer_start'      => 'inner_end',
          'seq_region_start' => 'seq_region_end',
          'inner_start'      => 'outer_end',

          'inner_end'        => 'outer_start',
          'seq_region_end'   => 'seq_region_start',
          'outer_end'        => 'inner_start',
        };
        my @order = qw/outer_start seq_region_start inner_start inner_end seq_region_end outer_end/;
        my $after_swap_mappings = {};
        foreach my $c (@order) {
          $final_mappings->{$id}->{$c} = $self->get_start($filtered_mappings->{$id}->{$swap_map->{$c}});
        }
        next;
      }
      foreach my $c (@new_coord_types) {
        my $start = $self->get_start($filtered_mappings->{$id}->{$c});
        $final_mappings->{$id}->{$c} = $start;
      }
    } else {
      $failed_mappings->{$id} = 'Incomplete mappings';
    }
  }

  my $file_filtered_mappings = $self->param('file_filtered_mappings');
  my $fh_filtered_mappings = FileHandle->new($file_filtered_mappings, 'w');

  foreach my $id (keys %$final_mappings) {
    my @values = ();
    my $seq_region_name = $lookup->{$id}->{seq_region_name};
    push @values, "seq_region_name=$seq_region_name";
    push @values, "structural_variation_feature_id=$id";
    foreach my $coord (keys %{$final_mappings->{$id}}) {
      my $start = $final_mappings->{$id}->{$coord};
      push @values, "$coord=$start";
    }
    print $fh_filtered_mappings join("\t", @values), "\n";
  }

  $fh_filtered_mappings->close();
}

sub join_svf_data {
  my $self = shift;

  my $file_init_feature = $self->param('file_init_feature');
  my $fh_init_feature = FileHandle->new($file_init_feature, 'r'); 

  my $feature_data = {};
  my $key = '';
  while (<$fh_init_feature>) {
    chomp;
    my $data = $self->read_data($_);
    $key = $data->{structural_variation_feature_id};
    $feature_data->{$key} = $data;
  }
  $fh_init_feature->close();

  # get new seq_region_ids
  my $seq_region_ids = {};
  my $cdba = $self->param('cdba_newasm');
  my $sa = $cdba->get_SliceAdaptor;
  my $slices = $sa->fetch_all('toplevel', undef, 1);
  foreach my $slice (@$slices) {
    $seq_region_ids->{$slice->seq_region_name} = $slice->get_seq_region_id;
  }
  # new map_weights
  my $file_filtered_mappings = $self->param('file_filtered_mappings');

  # join feature data with mapping data:
  my $file_load_features = $self->param('file_load_features');
  my $fh_load_features = FileHandle->new($file_load_features, 'w');   
  my $fh_mappings = FileHandle->new($file_filtered_mappings, 'r');
  my ($data, $variation_feature_id, $version, $variation_name);
  while (<$fh_mappings>) {
    chomp;
    my $filtered_data = $self->read_data($_); 
    my $seq_region_name = $filtered_data->{seq_region_name}; 
    my $seq_region_id = $seq_region_ids->{$seq_region_name};
    my $svf_id = $filtered_data->{structural_variation_feature_id}; 
    $data = $feature_data->{$svf_id};
    foreach my $coords (qw/outer_start seq_region_start inner_start inner_end seq_region_end outer_end/) {
      my $location = $data->{coord};
      $data->{seq_region_id} = $seq_region_id || '\N';
    }
    my $line = $self->print_complete_feature_line($data);
    print $fh_load_features $line, "\n"; 
  }

  $fh_mappings->close();
  $fh_load_features->close();
}



1;
