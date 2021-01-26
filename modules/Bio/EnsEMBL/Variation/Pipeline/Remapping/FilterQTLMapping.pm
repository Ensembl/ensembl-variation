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

package Bio::EnsEMBL::Variation::Pipeline::Remapping::FilterQTLMapping;

use strict;
use warnings;

use FileHandle;
use Bio::EnsEMBL::Registry;

use base ('Bio::EnsEMBL::Variation::Pipeline::Remapping::FilterMapping');

sub fetch_input {
  my $self = shift;
  $self->SUPER::fetch_input;
}

sub run {
  my $self = shift;
  $self->report_failed_qtl_mappings();
  $self->filter_qtl_mapping_results();
  $self->join_qtl_data();
  $self->SUPER::write_statistics();
}

sub write_output {
  my $self = shift;
}

sub report_failed_qtl_mappings {
  my $self = shift;
  $self->SUPER::report_failed_read_mappings;
}

sub filter_qtl_mapping_results {
  my $self = shift;
  $self->SUPER::filter_read_mapping_results('phenotype_feature_id');
}

sub join_qtl_data {
  my $self = shift;

  my $file_init_feature = $self->param('file_init_feature');
  my $fh_init_feature = FileHandle->new($file_init_feature, 'r'); 

  my $feature_data = {};
  while (<$fh_init_feature>) {
    chomp;
    my $data = $self->read_data($_);
    my $key = $data->{phenotype_feature_id};
    $feature_data->{$key} = $data;
  }
  $fh_init_feature->close();

  my $new_seq_region_ids = $self->get_new_seq_region_ids;

  # join feature data with mapping data:
  my $file_load_features = $self->param('file_load_features');
  my $fh_load_features = FileHandle->new($file_load_features, 'w');   
  my $file_filtered_mappings = $self->param('file_filtered_mappings');
  my $fh_mappings = FileHandle->new($file_filtered_mappings, 'r');
  my ($data, $variation_feature_id, $version, $variation_name);
  while (<$fh_mappings>) {
    chomp;
    my ($query_id, $seq_name, $start, $end, $strand, $score) = split("\t", $_);
    my $seq_region_id = $new_seq_region_ids->{$seq_name};
    $data = $feature_data->{$query_id};
    $data->{seq_region_id} = $seq_region_id;
    $data->{seq_region_start} = $start;
    $data->{seq_region_end} = $end;
    $data->{seq_region_strand} = $strand;
    my $line = $self->print_complete_feature_line($data);
    print $fh_load_features $line, "\n"; 
  }

  $fh_mappings->close();
  $fh_load_features->close();
}


1;
