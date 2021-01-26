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
 developers list at <https://lists.ensembl.org/mailman/listinfo/dev>.
 Questions may also be sent to the Ensembl help desk at
 <https://www.ensembl.org/Help/Contact>.
=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::Pipeline::Remapping::StructuralVariationFeatureQC;

use base ('Bio::EnsEMBL::Variation::Pipeline::Remapping::BaseRemapping');

use Bio::EnsEMBL::Registry;
use FileHandle;

sub fetch_input {
  my $self = shift;
  $self->SUPER::fetch_input;
}

sub run {
  my $self = shift;
  my $failed_sv_ids = $self->unmapped_structural_variation_features;
  my $prev_failed_sv_ids = $self->prev_failed_structural_variations; 
  $self->load_failed_sv($failed_sv_ids, $prev_failed_sv_ids);
  $self->rename_mapped_feature_table;
}

sub load_failed_sv {
  my $self = shift;
  my $failed_sv_ids = shift;
  my $prev_failed_sv_ids = shift;
  my $vdba_newasm = $self->get_newasm_variation_database_connection;
  my $dbc  = $vdba_newasm->dbc;
  my $dbh = $dbc->db_handle;

  if (!$self->table_exists('before_remapping_failed_structural_variation', $vdba_newasm) && $self->table_exists('failed_structural_variation', $vdba_newasm)) {
    $dbh->do("RENAME TABLE failed_structural_variation to before_remapping_failed_structural_variation;") or die $!;
    $dbh->do("CREATE TABLE failed_structural_variation LIKE before_remapping_failed_structural_variation;") or die $!;
  }
  if ($self->table_exists('before_remapping_failed_structural_variation', $vdba_newasm) && $self->table_exists('failed_structural_variation', $vdba_newasm)) {
    $dbh->do("TRUNCATE TABLE failed_structural_variation;") or die $!;
  }

  foreach my $sv_id (keys %$failed_sv_ids) {
    foreach my $failed_description_id (keys %{$failed_sv_ids->{$sv_id}}) {
      $dbh->do("INSERT INTO failed_structural_variation(structural_variation_id, failed_description_id) VALUES($sv_id, $failed_description_id);") or die $!;
    }
  }

  foreach my $sv_id (keys %$prev_failed_sv_ids) {
    foreach my $failed_description_id (keys %{$prev_failed_sv_ids->{$sv_id}}) {
      $dbh->do("INSERT INTO failed_structural_variation(structural_variation_id, failed_description_id) VALUES($sv_id, $failed_description_id);") or die $!;
    }
  }
}

sub prev_failed_structural_variations {
  my $self = shift;
  my $vdba = $self->get_oldasm_variation_database_connection;
  my $dbh = $vdba->dbc->db_handle;
  my $failed_variations = {};
  my $sth = $dbh->prepare(qq{
    SELECT structural_variation_id, failed_description_id FROM failed_structural_variation WHERE failed_description_id not in (17, 18);
  }, {mysql_use_result => 1});
  $sth->execute();
  while (my $row = $sth->fetchrow_arrayref) {
    $failed_variations->{$row->[0]}->{$row->[1]} = 1;
  }
  $sth->finish();
  return $failed_variations;
}

sub unmapped_structural_variation_features {
  my $self = shift;
  my $failed_sv_ids = {};

  my $dump_features_dir = $self->param('dump_features_dir');
  my $load_features_dir = $self->param('load_features_dir');

  # failed_descriptions TODO check that ids still match descriptions
  #  '17' => 'Variant can not be re-mapped to the current assembly',
  #  '18' => 'Supporting evidence can not be re-mapped to the current assembly',

  opendir (DIR, $load_features_dir) or die $!;
  while (my $load_file = readdir(DIR)) {
    next if ($load_file =~ /^\./);
    my $prev_mappings = _read_init_sv_info("$dump_features_dir/$load_file");
    my $new_mappings =  $self->_read_mapped_sv_info("$load_features_dir/$load_file");
    foreach my $sv_id (keys %$prev_mappings) {
      if (defined $new_mappings->{$sv_id}) {
        my ($has_mapping_prev, $has_evidence_prev) = @{_has_mapping_has_evidence($prev_mappings->{$sv_id})};
        my ($has_mapping_new, $has_evidence_new) = @{_has_mapping_has_evidence($new_mappings->{$sv_id})};
        if ($has_mapping_prev && !$has_mapping_new) {
          # lost mapping
          $failed_sv_ids->{$sv_id}->{17} = 1;
        }
        if ($has_evidence_prev && !$has_evidence_new) {
          # lost evidence
          $failed_sv_ids->{$sv_id}->{18} = 1;
        }
      } else {
        $failed_sv_ids->{$sv_id}->{17} = 1;
      }
    }
  }
  closedir (DIR);

  return $failed_sv_ids;
}

sub _has_mapping_has_evidence {
  my $mappings = shift;
  my $has_mapping = 0;
  my $has_evidence = 0;
  foreach my $svf_id (keys %$mappings) {
    if ($mappings->{$svf_id} eq 'is_mapping') {
      $has_mapping = 1;
    } else {
      $has_evidence = 1;
    }
  }
  return [$has_mapping, $has_evidence];  
}

sub _read_mapped_sv_info {
  my $self = shift;
  my $file = shift;

  my $vdba_newasm = $self->get_newasm_variation_database_connection;
  my @columns = @{$self->get_sorted_column_names($vdba_newasm, 'structural_variation_feature')};

  my $sv_info = {}; 
  my $fh = FileHandle->new($file, 'r');

  while (<$fh>) {
    chomp;
    my @values = split("\t", $_);
    my %data = map {$columns[$_] => $values[$_]} (0..(scalar @columns - 1));
    my $sv_id = $data{'structural_variation_id'};
    my $svf_id = $data{'structural_variation_feature_id'};
    $sv_info->{$sv_id}->{$svf_id} = ($data{'is_evidence'}) ? 'is_evidence' : 'is_mapping';
  }

  return $sv_info;
}
sub _read_init_sv_info {
  my $file = shift;
  my $sv_info = {}; 
  my $fh = FileHandle->new($file, 'r');

  while (<$fh>) {
    chomp;
    my @key_values = split("\t", $_);
    my $mapping = {};
    foreach my $key_value (@key_values) {
      my ($column_name, $value) = split('=', $key_value, 2);
      $mapping->{$column_name} = $value;
    }
    my $sv_id = $mapping->{'structural_variation_id'};
    my $svf_id = $mapping->{'structural_variation_feature_id'};
    $sv_info->{$sv_id}->{$svf_id} = ($mapping->{'is_evidence'}) ? 'is_evidence' : 'is_mapping';
  }
  return $sv_info;
}

1;
