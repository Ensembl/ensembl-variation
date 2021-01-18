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

package Bio::EnsEMBL::Variation::Pipeline::Remapping::InitVariationFeatureQC;
use base ('Bio::EnsEMBL::Variation::Pipeline::Remapping::BaseRemapping');

use Bio::DB::Fasta;
use Bio::Perl;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp expand);
use FileHandle;

sub fetch_input {
  my $self = shift;
  $self->SUPER::fetch_input;
}

sub run {
  my $self = shift;
  my $vdba = $self->get_newasm_variation_database_connection;
  my $dbh = $vdba->dbc->db_handle();

  my $qc_mapped_features_dir = $self->param('qc_mapped_features_dir');
  my $file_count = 1;
  my $entries_per_file = $self->param('entries_per_file');
  my $count_entries = 0;
  my $fh = FileHandle->new("$qc_mapped_features_dir/$file_count.txt", 'w');

  my @column_names = qw/variation_feature_id seq_region_id seq_region_name seq_region_start seq_region_end seq_region_strand allele_string variation_id variation_name map_weight alignment_quality/;
  my $column_concat = join(',', @column_names);

  my $feature_table = $self->param('feature_table') . '_mapping_results';

  my $sth = $dbh->prepare(qq{
    SELECT vf.variation_feature_id, vf.seq_region_id, sr.name, vf.seq_region_start, vf.seq_region_end, vf.seq_region_strand, vf.allele_string, vf.variation_id, vf.variation_name, vf.map_weight, vf.alignment_quality
    FROM $feature_table vf, seq_region sr
    WHERE sr.seq_region_id = vf.seq_region_id;
  }, {mysql_use_result => 1});

  $sth->execute();
  while (my $row = $sth->fetchrow_arrayref) {
    my @values = map { defined $_ ? $_ : '\N'} @$row;
    my @pairs = ();
    for my $i (0..$#column_names) {
      push @pairs, "$column_names[$i]=$values[$i]";
    }
    if ($count_entries >= $entries_per_file) {
      $fh->close();
      $file_count++;
      $fh = FileHandle->new("$qc_mapped_features_dir/$file_count.txt", 'w');
      $count_entries = 0;
    }
    $count_entries++;
    print $fh join("\t", @pairs), "\n";
  }

  $sth->finish();
  $fh->close();

  $self->param('file_count', $file_count);

  my $fasta_db_dir = $self->param('new_assembly_dir');
  $self->test_bioperl_version;
  my $fasta_db = Bio::DB::Fasta->new($fasta_db_dir, -reindex => 1);
  $self->param('fasta_db', $fasta_db_dir);

}

sub write_output {
  my $self = shift;
  my @jobs = ();
  my $file_count = $self->param('file_count');
  my $i = 1;
  while ($i <= $file_count) {
    push @jobs, {
      'file_number' => $i,
      'fasta_db' => $self->param('fasta_db'),
    };
    $i++;
  } 
  $self->dataflow_output_id(\@jobs, 2)
}

1;
