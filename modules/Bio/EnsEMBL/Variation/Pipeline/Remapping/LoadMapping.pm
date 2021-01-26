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
package Bio::EnsEMBL::Variation::Pipeline::Remapping::LoadMapping;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use ImportUtils qw(load);
use base ('Bio::EnsEMBL::Variation::Pipeline::Remapping::BaseRemapping');

sub fetch_input {
  my $self = shift;
  $self->SUPER::fetch_input;
}

sub run {
  my $self = shift;
  if ($self->param('mode') eq 'remap_read_coverage') {
    $self->load_read_coverage();
  } elsif ($self->param('mode') eq 'remap_QTL' || $self->param('mode') eq 'remap_svf') {
    $self->load_mapping_results();
  } else {
    $self->load_features();
  }
}

sub write_output {
  my $self = shift;
}

sub load_read_coverage {
  my $self = shift;

  my $load_features_dir = $self->param('load_features_dir');
  my $vdba = $self->get_newasm_variation_database_connection;
  my $dbc = $vdba->dbc;

  my $feature_table = $self->param('feature_table');
  my $result_table = "$feature_table\_mapping_results";

  $dbc->do(qq{ DROP TABLE IF EXISTS $result_table});
  $dbc->do(qq{ CREATE TABLE $result_table like $feature_table });
  $dbc->do(qq{ ALTER TABLE $result_table DISABLE KEYS});

  my @column_names = @{$self->get_sorted_column_names($vdba, $result_table)};
  my $column_names_concat = join(',', @column_names);

  opendir(IND_DIR, $load_features_dir) or die $!;
  while (my $individual_dir = readdir(IND_DIR))  {
    next if ($individual_dir =~ /^\./);
    my $dir = "$load_features_dir/$individual_dir";
    opendir(DIR, $dir) or die $!;
    while (my $file = readdir(DIR)) {
      if ($file =~ /^(.+)\.txt$/) {
        my $file_number = $1;
        $self->run_cmd("cp $dir/$file_number.txt $dir/load_$file_number.txt");
        my $TMP_DIR = $dir;
        my $tmp_file = "load_$file_number.txt";
        $ImportUtils::TMP_DIR = $TMP_DIR;
        $ImportUtils::TMP_FILE = $tmp_file;
        load($dbc, ($result_table, $column_names_concat));
      }
    }
    closedir(DIR);
  }   
  closedir(IND_DIR);

  $dbc->do(qq{ ALTER TABLE $result_table ENABLE KEYS});
}

sub load_mapping_results {
  my $self = shift;

  my $vdba = $self->get_newasm_variation_database_connection;
  my $dbc  = $vdba->dbc;

  my $load_features_dir = $self->param('load_features_dir');
  my $feature_table = $self->param('feature_table');
  my $result_table = $self->param('feature_table_mapping_results');

  $dbc->do(qq{ DROP TABLE IF EXISTS $result_table});
  $dbc->do(qq{ CREATE TABLE $result_table like $feature_table });
  $dbc->do(qq{ ALTER TABLE $result_table DISABLE KEYS});

  my @column_names = @{$self->get_sorted_column_names($vdba, $result_table)};
  my $column_names_concat = join(',', @column_names);
  opendir(DIR, $load_features_dir) or die $!;
  while (my $file = readdir(DIR)) {
    if ($file =~ /^(.+)\.txt$/) {
      my $file_number = $1;
      $self->run_cmd("cp $load_features_dir/$file_number.txt $load_features_dir/load_$file_number.txt");
      my $TMP_DIR = $load_features_dir;
      my $tmp_file = "load_$file_number.txt";
      $ImportUtils::TMP_DIR = $TMP_DIR;
      $ImportUtils::TMP_FILE = $tmp_file;
      load($dbc, ($result_table, $column_names_concat));
    }
  }
  closedir(DIR);
  $dbc->do(qq{ ALTER TABLE $result_table ENABLE KEYS});
}

sub load_features {
  my $self = shift;

  my $vdba = $self->get_newasm_variation_database_connection;
  my $dbc  = $vdba->dbc;

  my $load_features_dir = $self->param('load_features_dir');
  my $feature_table = $self->param('feature_table');
  my $result_table = $self->param('feature_table_mapping_results');

  $dbc->do(qq{ DROP TABLE IF EXISTS $result_table});
  $dbc->do(qq{ CREATE TABLE $result_table like $feature_table });
  $dbc->do(qq{ ALTER TABLE $result_table ADD variation_feature_id_old int(10) unsigned});
  $dbc->do(qq{ ALTER TABLE $result_table DISABLE KEYS});

  my @column_names = @{$self->get_sorted_column_names($vdba, $result_table, {variation_feature_id => 1})};
  my $column_names_concat = join(',', @column_names);

  opendir(DIR, $load_features_dir) or die $!;
  while (my $file = readdir(DIR)) {
    if ($file =~ /^(.+)\.txt$/) {
      my $file_number = $1;
      $self->run_cmd("cp $load_features_dir/$file_number.txt $load_features_dir/load_$file_number.txt");
      my $TMP_DIR = $load_features_dir;
      my $tmp_file = "load_$file_number.txt";
      $ImportUtils::TMP_DIR = $TMP_DIR;
      $ImportUtils::TMP_FILE = $tmp_file;
      load($dbc, ($result_table, $column_names_concat));
    }
  }
  closedir(DIR);

  $dbc->do(qq{ ALTER TABLE $result_table ENABLE KEYS});
}

1;
