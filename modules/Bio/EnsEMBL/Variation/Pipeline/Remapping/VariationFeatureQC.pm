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

use strict;
use warnings;

package Bio::EnsEMBL::Variation::Pipeline::Remapping::VariationFeatureQC;

use base ('Bio::EnsEMBL::Variation::Pipeline::Remapping::BaseRemapping');

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp expand);
use Bio::EnsEMBL::Variation::Utils::RemappingUtils qw(qc_mapped_vf);
use FileHandle;

sub fetch_input {
  my $self = shift;
  $self->SUPER::fetch_input;
}

sub run {
  my $self = shift;
  my $file_number = $self->param('file_number');
  my $fasta_db = $self->param('fasta_db');
  my $vdba = $self->get_newasm_variation_database_connection;
  my $vdba_oldasm = $self->get_oldasm_variation_database_connection;
  my $qc_mapped_features_dir = $self->param('qc_mapped_features_dir');
  my $qc_update_features_dir = $self->param('qc_update_features_dir');
  my $qc_failure_reasons_dir = $self->param('qc_failure_reasons_dir');

  my $feature_table = $self->param('feature_table') . '_mapping_results';
  my $config = {
    fasta_db => $fasta_db,
    mapped_features_file => "$qc_mapped_features_dir/$file_number.txt",
    update_features_file => "$qc_update_features_dir/$file_number.txt",
    failure_reasons_file => "$qc_failure_reasons_dir/$file_number.txt",
    feature_table => $feature_table,
    vdba => $vdba,
    vdba_oldasm => $vdba_oldasm,
  };
  qc_mapped_vf($config);

  # Check that there are no duplicated entries in variation_feature after remapping

  my $dbh = $vdba->dbc->db_handle();
  my $sth = $dbh->prepare(qq{
    SELECT count(vf1.variation_feature_id)
    FROM $feature_table vf1
    JOIN $feature_table vf2
    ON (vf2.variation_id = vf1.variation_id
      AND vf2.variation_feature_id > vf1.variation_feature_id
      AND vf2.seq_region_id = vf1.seq_region_id
      AND vf2.seq_region_start = vf1.seq_region_start
      AND vf2.seq_region_end = vf1.seq_region_end
      AND vf2.seq_region_strand = vf1.seq_region_strand
    );
  }, {mysql_use_result => 1});

  $sth->execute() or die $sth->errstr;
  my $row = $sth->fetchrow_arrayref;
  $sth->finish();
  if ($row->[0] > 0) {
    die "$feature_table contains duplicates after remapping.\n";
  }
}

1;
