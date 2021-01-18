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
package Bio::EnsEMBL::Variation::Pipeline::Remapping::BaseRemapping;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Hive::Process);

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Registry;
use Bio::Perl;

sub get_newasm_variation_database_connection {
  my $self = shift;
  return $self->_get_database_connection('variation', 'registry_file_newasm');
}

sub get_oldasm_variation_database_connection {
  my $self = shift;
  return $self->_get_database_connection('variation', 'registry_file_oldasm');
}

sub get_newasm_core_database_connection {
  my $self = shift;
  return $self->_get_database_connection('core', 'registry_file_newasm');
}

sub get_oldasm_core_database_connection {
  my $self = shift;
  return $self->_get_database_connection('core', 'registry_file_oldasm');
}

sub _get_database_connection {
  my $self = shift;
  my $type = shift;
  my $registry_file = shift;
  my $registry = 'Bio::EnsEMBL::Registry';
  $registry->load_all($self->param($registry_file));
  return  $registry->get_DBAdaptor($self->param('species'), $type);
}

sub test_bioperl_version {
  my $self = shift;
  my $bioperl_version = Bio::Perl->VERSION;
  die "At least Bio::Perl 1.006924 required" if ($bioperl_version < 1.006924);
}

sub read_line {
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

sub count_files {
  my $self = shift;
  my $dir = shift;
  my $file_type = shift;
  opendir(my $dh, $dir) or die "Not a directory $dir";
  my $count = scalar(grep { $_ =~ m/$file_type$/ } readdir($dh));
  closedir($dh);
  return $count;
}

sub run_cmd {
  my $self = shift;
  my $cmd = shift;
  if (my $return_value = system($cmd)) {
    $return_value >>= 8;
    die "system($cmd) failed: $return_value";
  }
}

sub get_individual_name {
  my $self = shift;
  my $individual_id = shift;
  my $vdba = $self->param('vdba_oldasm');
  my $dbname = $vdba->dbc->dbname();
  my $dbh = $vdba->dbc->db_handle();

  my $column = 'sample_id';
  my $table = 'sample';
  my $column_names = $self->get_column_names($dbh, $dbname, 'individual');
  if (grep /individual_id/, @$column_names) {
    $column = 'individual_id';
    $table = 'individual';
  }

  my $sth = $dbh->prepare(qq{SELECT name FROM $table where $column=$individual_id;});
  my @names = ();
  $sth->execute();
  while (my $row = $sth->fetchrow_arrayref) {
    push @names, $row->[0];
  }
  $sth->finish();
  die ("Wrong number of names for individiual_id: $individual_id in DB: $dbname: " . scalar @names) if ((scalar @names) != 1);
  return $names[0];
}

sub get_column_names {
  my ($self, $dbh, $dbname, $table_name) = @_;
  my $query = qq{SHOW columns FROM $dbname.$table_name};
  my $column_names_info = $self->run_query($dbh, $query);
  my @column_names = ();
  foreach my $column_name_info (@$column_names_info) {
    my ($name, $info) = split(',', $column_name_info, 2);
    push @column_names, $name;
  }
  return \@column_names;
}

sub get_seq_region_ids {
  my ($self, $cdba) = @_;
  my $sa = $cdba->get_SliceAdaptor;
  # don't include asm exceptions
  my $slices = $sa->fetch_all('toplevel', undef, 0, 1);
  my $seq_region_ids = {};
  foreach my $slice (@$slices) {
    my $seq_region_name = $slice->seq_region_name;
    next if ($seq_region_name =~ /PATCH/);
    my $seq_region_id = $slice->get_seq_region_id;
    $seq_region_ids->{$seq_region_id} = $seq_region_name;
  }
  return $seq_region_ids;
}

sub get_table_names {
  my ($self, $dbh, $dbname) = @_;
  my $query = qq{
    SHOW TABLES FROM $dbname;
  };
  my $table_names = $self->run_query($dbh, $query);
  return $table_names;
}

sub run_query {
  my ($self, $dbh, $query) = @_;
  my $sth = $dbh->prepare($query);
  $sth->execute();
  my @results = ();
  while (my $row = $sth->fetchrow_arrayref) {
    my @values = map { defined $_ ? $_ : '\N' } @$row;
    push @results, join(',', @values);
  }
  $sth->finish();
  return \@results;
}

sub get_sorted_column_names {
  my $self = shift;
  my $dba = shift;
  my $table = shift;
  my $exclude = shift;
  my $dbc  = $dba->dbc;
  my $dbname = $dba->dbc->dbname();
  my $dbh = $dbc->db_handle;
  my $sth = $dbh->prepare(qq{
      SELECT COLUMN_NAME FROM INFORMATION_SCHEMA.COLUMNS
      WHERE TABLE_SCHEMA = '$dbname'
      AND TABLE_NAME = '$table';
      });
  $sth->execute();

  my @column_names = ();
  while (my @name = $sth->fetchrow_array) {
    if (!$exclude->{$name[0]}) {
      push @column_names, $name[0];
    }
  }
  $sth->finish();
  @column_names = sort @column_names;
  return \@column_names;
}

sub table_exists {
  my $self = shift;
  my $table_name = shift;
  my $vdba = shift;
  my $dbc  = $vdba->dbc;
  my $dbh = $dbc->db_handle;
  my $sth = $dbh->prepare(qq{
    SHOW TABLES LIKE ?;
  });
  $sth->execute($table_name);
  my @rows = @{$sth->fetchall_arrayref};
  $sth->finish;
  return scalar @rows;
}

sub rename_mapped_feature_table {
  my $self = shift;
  my $feature_table = $self->param('feature_table');
  my $result_table = $self->param('feature_table_mapping_results');

  my $vdba_newasm = $self->get_newasm_variation_database_connection;
  my $dbh  = $vdba_newasm->dbc->db_handle;

  # only create before_remapping table once
  if (!$self->table_exists("before_remapping_$feature_table", $vdba_newasm) && $self->table_exists($feature_table, $vdba_newasm)) {
    $dbh->do("RENAME TABLE $feature_table TO before_remapping_$feature_table;") or die $!;
  }

  if ($self->table_exists($result_table, $vdba_newasm)) {
    # drop feature table if it exists
    # if it exists at this point it means that before_remapping_feature_table has been created and this is another remapping run
    $dbh->do("DROP TABLE IF EXISTS $feature_table;") or die $!;
    $dbh->do("RENAME TABLE $result_table TO $feature_table;") or die $!;
  }
}

1;
