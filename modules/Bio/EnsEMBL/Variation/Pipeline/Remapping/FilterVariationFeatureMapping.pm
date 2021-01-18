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

package Bio::EnsEMBL::Variation::Pipeline::Remapping::FilterVariationFeatureMapping;

use strict;
use warnings;

use FileHandle;
use Bio::EnsEMBL::Registry;

use base ('Bio::EnsEMBL::Variation::Pipeline::Remapping::FilterMapping');
use Bio::EnsEMBL::Variation::Utils::RemappingUtils qw(filter_vf_mapping_first_round filter_vf_mapping_second_round);

# failure_reasons
my $NO_MAPPING = 'no mapping';
my $TOO_MANY_MAPPINGS = 'map weight > 5';
my $POOR_SCORE = 'does not pass alignment score';

sub fetch_input {
  my $self = shift;
  $self->SUPER::fetch_input;
}

sub run {
  my $self = shift;
  $self->report_failed_mappings();
  $self->filter_mapping_results();
  $self->join_feature_data();
  $self->SUPER::write_statistics();
}

sub write_output {
  my $self = shift;
}

sub report_failed_mappings {
  my $self = shift;
  my $file_mappings = $self->param('file_mappings');
  my $file_failed_mappings = $self->param('file_failed_mappings');

  my $mapped = {};
  my $unmapped = {};
  my $fh_mappings = FileHandle->new($file_mappings, 'r');
  while (<$fh_mappings>) {
    chomp;
    my ($old_seq_info, $new_seq_info, $query_name, $map_weight, $cigar, $relative_algn_score, $clipped_info) = split("\t", $_);
    $mapped->{$query_name} = 1;
  }
  $fh_mappings->close();
  
  my $fh_failed_mappings = FileHandle->new($file_failed_mappings, 'r');
  while (<$fh_failed_mappings>) {
    chomp;
    my ($indel, $old_seq_info, $new_seq_info, $query_name, $map_weight, $cigar, $relative_algn_score, $clipped_info) = split("\t", $_);
    unless ($mapped->{$query_name}) {
      $unmapped->{$query_name} = 1;
    }
  }
  $fh_failed_mappings->close();

  my $file_failed_filtered_mappings = $self->param('file_failed_filtered_mappings');
  my $fh_failed_filtered_mappings = FileHandle->new($file_failed_filtered_mappings, 'w');
  foreach my $query_name (keys %$unmapped) {
    print $fh_failed_filtered_mappings "$query_name\t$NO_MAPPING\n";
  }
  $fh_failed_filtered_mappings->close();

  my $count_mapped = scalar keys %$mapped;
  my $count_unmapped = scalar keys %$unmapped;

  $self->param('pre_count_mapped', $count_mapped);
  $self->param('pre_count_unmapped', $count_unmapped);

  # get input id counts
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

sub filter_mapping_results {
  my $self = shift;
  my $algn_score_threshold = $self->param('algn_score_threshold');
  my $use_prior_info = $self->param('use_prior_for_filtering');
  my $map_to_chrom_only = $self->param('map_to_chrom_only');
  my $max_map_weight = $self->param('max_map_weight');
  my $chroms = $self->get_seq_region_names('chromosome') if ($map_to_chrom_only);
  my $count_chroms = keys %$chroms;

  my $file_mappings = $self->param('file_mappings');
  my $file_filtered_mappings = $self->param('file_filtered_mappings');
  my $file_failed_filtered_mappings = $self->param('file_failed_filtered_mappings');

  my $stats_config = {};  
  my $config = {};
  foreach (qw/algn_score_threshold use_prior_for_filtering map_to_chrom_only max_map_weight file_mappings file_filtered_mappings file_failed_filtered_mappings/) {
    $config->{$_} = $self->param($_) || undef;
  }

  my $filter_results_first_round =  filter_vf_mapping_first_round($config, $stats_config, $chroms);

  my $multi_map_all = $filter_results_first_round->{multi_map_all};
  my $multi_map_same_chrom = $filter_results_first_round->{multi_map_same_chrom};

  my $count_multi_map_all        = scalar keys %$multi_map_all;
  my $count_multi_map_same_chrom = scalar keys %$multi_map_same_chrom;

  my $stats_failed_non_chrom = $count_multi_map_all - $count_multi_map_same_chrom;

  my $multi_map_working = ($use_prior_info) ? $multi_map_same_chrom : $multi_map_all;

  filter_vf_mapping_second_round($config, $stats_config, $multi_map_working);

  $self->param('stats_failed_non_chrom', $stats_failed_non_chrom);
  foreach my $stats_value (qw/stats_unique_map stats_unique_map_after_filter stats_multi_map stats_multi_map_after_filter stats_failed_poor_score stats_failed_after_filter stats_exceeds_max_map_weight/) {
    $self->param($stats_value, $stats_config->{$stats_value});
  }
}

sub join_feature_data {
  my $self = shift;

  my $file_init_feature = $self->param('file_init_feature');
  my $fh_init_feature = FileHandle->new($file_init_feature, 'r'); 

  my $feature_data = {};
  while (<$fh_init_feature>) {
    chomp;
    my $data = $self->read_data($_);
    $feature_data->{$data->{variation_feature_id}} = $data;
  }
  $fh_init_feature->close();

  # get new seq_region_ids
  my $seq_region_ids = {};
  my $cdba = $self->get_newasm_core_database_connection;
  my $sa = $cdba->get_SliceAdaptor;
  my $slices = $sa->fetch_all('toplevel', undef, 1);
  foreach my $slice (@$slices) {
    $seq_region_ids->{$slice->seq_region_name} = $slice->get_seq_region_id;
  }
  # new map_weights
  my $file_filtered_mappings = $self->param('file_filtered_mappings');
  my $fh_mappings = FileHandle->new($file_filtered_mappings, 'r');
  my $map_weights = {};
  while (<$fh_mappings>) {
    chomp;
    my ($query_name, $seq_name, $start, $end, $strand, $score) = split("\t", $_);
    $map_weights->{$query_name}++;
  }
  $fh_mappings->close();

  # join feature data with mapping data:

  my $file_load_features = $self->param('file_load_features');
  my $fh_load_features = FileHandle->new($file_load_features, 'w');   
  $fh_mappings = FileHandle->new($file_filtered_mappings, 'r');
  my ($data, $variation_feature_id, $version, $variation_name);
  while (<$fh_mappings>) {
    chomp;
    my ($query_name, $seq_name, $start, $end, $strand, $score) = split("\t", $_);
    # query_name: 156358-150-1-150-11:5502587:5502587:1:T/C:rs202026261:dbSNP:SNV
    my @query_name_components = split('-', $query_name, 2);
    $variation_feature_id = $query_name_components[0];
    my $seq_region_id = $seq_region_ids->{$seq_name};
    my $map_weight = $map_weights->{$query_name};
    $data = $feature_data->{$variation_feature_id};
    $data->{seq_region_id} = $seq_region_id;
    $data->{seq_region_start} = $start;
    $data->{seq_region_end} = $end;
    $data->{seq_region_strand} = $strand;
    $data->{alignment_quality} = $score;
    $data->{map_weight} = $map_weight;

    my $line = $self->print_complete_feature_line($data);

    print $fh_load_features $line, "\n"; 
  }

  $fh_mappings->close();
  $fh_load_features->close();
}

1;
