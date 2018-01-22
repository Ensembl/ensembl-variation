#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.




=head1 CONTACT

Please email comments or questions to the public Ensembl
developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

Questions may also be sent to the Ensembl help desk at
<helpdesk.org>.

=cut
package Bio::EnsEMBL::Variation::Pipeline::Remapping::LoadMapping;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use ImportUtils qw(load);
use base ('Bio::EnsEMBL::Hive::Process');

sub fetch_input {
  my $self = shift;
  my $registry = 'Bio::EnsEMBL::Registry';
  $registry->load_all($self->param('registry_file_newasm'));
  my $vdba = $registry->get_DBAdaptor($self->param('species'), 'variation');
  $self->param('vdba', $vdba);
}

sub run {
  my $self = shift;
  if ($self->param('mode') eq 'remap_read_coverage') {
    $self->load_read_coverage();
  } elsif ($self->param('mode') eq 'remap_qtls') {
    $self->load_qtl();
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
  my $vdba = $self->param('vdba');
  my $dbc = $vdba->dbc;

  my $dbname = $vdba->dbc->dbname();
#  my $read_coverage_table = 'read_coverage';
#  my $result_table = "$read_coverage_table\_mapping_results";
  my $feature_table = $self->param('feature_table');
  my $result_table = "$feature_table\_mapping_results";


  $dbc->do(qq{ DROP TABLE IF EXISTS $result_table});
  $dbc->do(qq{ CREATE TABLE $result_table like $feature_table });
  $dbc->do(qq{ ALTER TABLE $result_table DISABLE KEYS});

  my $dbh = $dbc->db_handle;
  my $sth = $dbh->prepare(qq{
      SELECT COLUMN_NAME FROM INFORMATION_SCHEMA.COLUMNS
      WHERE TABLE_SCHEMA = '$dbname'
      AND TABLE_NAME = '$result_table';
      });
  $sth->execute();

# QC that all necessary columns are there: e.g. seq_region_id, ...
  my @column_names = ();
  while (my @name = $sth->fetchrow_array) {
    push @column_names, $name[0];
  }
  $sth->finish();
  @column_names = sort @column_names;

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
  $dbh->disconnect;
}

sub load_qtl {
  my $self = shift;

  my $vdba = $self->param('vdba');
  my $dbc  = $vdba->dbc;

  my $load_features_dir = $self->param('load_features_dir');

  my $dbname = $vdba->dbc->dbname();
  my $feature_table = $self->param('feature_table');
  my $result_table = "$feature_table\_mapping_results";
  $dbc->do(qq{ DROP TABLE IF EXISTS $result_table});
  $dbc->do(qq{ CREATE TABLE $result_table like $feature_table });
  $dbc->do(qq{ ALTER TABLE $result_table DISABLE KEYS});

  my $dbh = $dbc->db_handle;
  my $sth = $dbh->prepare(qq{
      SELECT COLUMN_NAME FROM INFORMATION_SCHEMA.COLUMNS
      WHERE TABLE_SCHEMA = '$dbname'
      AND TABLE_NAME = '$result_table';
      });
  $sth->execute();

# QC that all necessary columns are there: e.g. seq_region_id, ...
  my @column_names = ();
  while (my @name = $sth->fetchrow_array) {
    push @column_names, $name[0];
  }
  $sth->finish();
  @column_names = sort @column_names;

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
  $dbh->disconnect;
}
sub load_features {
  my $self = shift;

  my $vdba = $self->param('vdba');
  my $dbc  = $vdba->dbc;

  my $load_features_dir = $self->param('load_features_dir');

  my $dbname = $vdba->dbc->dbname();
  my $feature_table = $self->param('feature_table');
  my $result_table = "$feature_table\_mapping_results";
  if ($self->param('mode') eq 'remap_multi_map') {
    $result_table = "$feature_table\_multi_map_results";
  }
  if ($self->param('mode') eq 'remap_alt_loci') {
    $result_table = "$feature_table\_map_alt_loci_results";
  }
  if ($self->param('mode') eq 'remap_post_projection') {
    $result_table = "$feature_table\_post_projection";
  }
  $dbc->do(qq{ DROP TABLE IF EXISTS $result_table});
  $dbc->do(qq{ CREATE TABLE $result_table like $feature_table });
  $dbc->do(qq{ ALTER TABLE $result_table ADD variation_feature_id_old int(10) unsigned});
  $dbc->do(qq{ ALTER TABLE $result_table DISABLE KEYS});

  my $dbh = $dbc->db_handle;
  my $sth = $dbh->prepare(qq{
      SELECT COLUMN_NAME FROM INFORMATION_SCHEMA.COLUMNS
      WHERE TABLE_SCHEMA = '$dbname'
      AND TABLE_NAME = '$result_table';
      });
  $sth->execute();

# QC that all necessary columns are there: e.g. seq_region_id, ...
  my @column_names = ();
  while (my @name = $sth->fetchrow_array) {
    unless ($name[0] =~ /^variation_feature_id$/) {
      push @column_names, $name[0];
    }
  }
  $sth->finish();
  @column_names = sort @column_names;

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

  if ($self->param('mode') eq 'remap_post_projection') {
    my $feature_table_projection = $self->param('feature_table_projection');
    my $sth = $dbh->prepare(qq{
        SELECT COLUMN_NAME FROM INFORMATION_SCHEMA.COLUMNS
        WHERE TABLE_SCHEMA = '$dbname'
        AND TABLE_NAME = '$feature_table_projection';
        });
    $sth->execute();

# QC that all necessary columns are there: e.g. seq_region_id, ...
    my @column_names = ();
    while (my @name = $sth->fetchrow_array) {
      unless ($name[0] =~ /^variation_feature_id$/) {
        push @column_names, $name[0];
      }
    }
    $sth->finish();
    $column_names_concat = join(',' , @column_names);
    $dbc->do(qq{ INSERT INTO $result_table($column_names_concat) SELECT $column_names_concat FROM $feature_table_projection;});
  }
  $dbc->do(qq{ ALTER TABLE $result_table ENABLE KEYS});
  $dbh->disconnect;
}

sub run_cmd {
  my ($self, $cmd) = @_;
  if (my $return_value = system($cmd)) {
    $return_value >>= 8;
    die "system($cmd) failed: $return_value";
  }
}


1;
