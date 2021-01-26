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

package Bio::EnsEMBL::Variation::Pipeline::Remapping::FilterStructuralVariationFeatureMapping;

use strict;
use warnings;

use FileHandle;
use Bio::DB::Fasta;
use Bio::EnsEMBL::Registry;

use base ('Bio::EnsEMBL::Variation::Pipeline::Remapping::FilterMapping');
use Bio::EnsEMBL::Variation::Utils::RemappingUtils qw(filter_svf_mapping);

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

  my $file_mappings = $self->param('file_mappings');
  my $file_number = $self->param('file_number');
  my $dump_features_dir = $self->param('dump_features_dir'); 
  my $file_filtered_mappings = $self->param('file_filtered_mappings');
  my $file_failed_filtered_mappings = $self->param('file_failed_filtered_mappings');

  my $config = {
    file_prev_mappings => "$dump_features_dir/lookup_$file_number.txt",
    file_mappings => $file_mappings,
    file_filtered_mappings => $file_filtered_mappings,
    file_failed_filtered_mappings => "$file_failed_filtered_mappings",
  };
  filter_svf_mapping($config);
}

sub join_svf_data {
  my $self = shift;

  my $file_init_feature = $self->param('file_init_feature');
  my $fh_init_feature = FileHandle->new($file_init_feature, 'r'); 
  my $init_feature_data = {};
  while (<$fh_init_feature>) {
    chomp;
    my $data = $self->read_data($_);
    $init_feature_data->{$data->{structural_variation_feature_id}} = $data;
  }
  $fh_init_feature->close();

  my $new_seq_region_ids = $self->get_new_seq_region_ids;

  # join feature data with mapping data:
  my $file_load_features = $self->param('file_load_features');
  my $fh_load_features = FileHandle->new($file_load_features, 'w');   

  my $file_filtered_mappings = $self->param('file_filtered_mappings');
  my $fh_mappings = FileHandle->new($file_filtered_mappings, 'r');
  while (<$fh_mappings>) {
    chomp;
    my $filtered_data = $self->read_data($_); 
    my $seq_region_name = $filtered_data->{seq_region_name}; 
    my $seq_region_id = $new_seq_region_ids->{$seq_region_name};
    die "No new seq_region_id for $seq_region_name" unless($seq_region_id);
    my $svf_id = $filtered_data->{structural_variation_feature_id}; 
    my $data = $init_feature_data->{$svf_id};
    $data->{seq_region_id} = $seq_region_id;
    foreach my $coord (qw/outer_start seq_region_start inner_start inner_end seq_region_end outer_end seq_region_strand/) {
      $data->{$coord} = $filtered_data->{$coord} || '\N';
    }
    my $line = $self->print_complete_feature_line($data);
    print $fh_load_features $line, "\n"; 
  }

  $fh_mappings->close();
  $fh_load_features->close();
}

sub print_complete_feature_line {
  my $self = shift;
  my $data = shift;
  my @output = ();
  foreach my $column_name (sort keys %$data) {
    unless ($column_name eq 'seq_region_name') {
      push @output, $data->{$column_name};
    }
  }
  my $line = join("\t", @output);
  return $line;
}

1;
