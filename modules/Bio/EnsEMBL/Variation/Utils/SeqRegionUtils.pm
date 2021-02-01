=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2021] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     https://www.apache.org/licenses/LICENSE-2.0

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


=head1 SeqRegionUtils

This module is a general package based on the perl script update_seq_region_ids.pl

=cut

package Bio::EnsEMBL::Variation::Utils::SeqRegionUtils;

use strict;
use warnings;

use base qw(Exporter);

our @EXPORT_OK = qw(update_seq_region_ids);


=head2 update_seq_region_ids

  Arg[1]      : core DB adaptor
  Arg[2]      : variation DB adaptor
  Arg[3]      : (optional) dry_run: only list the sql updates that would be applied; default 0
  Example     :  use Bio::EnsEMBL::Variation::Pipeline::Utils::SeqRegionUtils qw(update_seq_region_ids);
                update_seq_region_ids($core_dba, $variation_dba, 0)
  Description : checks that seq_region name and seq_region_id are in sync
                with the core db and applies any needed updates.
  Exceptions  : None
  Caller      : Phenotype import pipeline for human import

=cut

sub update_seq_region_ids {
  my $core_dba = shift;
  my $variation_dba = shift;
  my $dry_run = shift // 0;

  my $dbname = $core_dba->dbc->dbname;
  my $core_species = $core_dba->species;
  my $variation_species = $variation_dba->species;
  return if ($variation_species ne $core_species);

  my $dbh = $core_dba->dbc->db_handle;
  my $max_mapping_set_id = get_max_mapping_set_id($dbh, $dbname);
  return if $max_mapping_set_id == 0;

  # check the mapping set contains applicable records
  my $sthCheck = $dbh->prepare("SELECT internal_schema_build, external_schema_build FROM mapping_set WHERE mapping_set_id=$max_mapping_set_id;");
  $sthCheck->execute();
  my @row = $sthCheck->fetchrow_array;
  my @name_parts = split('_', $dbname); # dbname is eg. bos_taurus_variation_104_12
  my $current_db_schema_version = $name_parts[-2] . "_". $name_parts[-1]; #$current_db_schema_version will be 104_12
  return if $row[0] ne $current_db_schema_version;

  my $id_mapping = {};
  my $sth = $dbh->prepare("SELECT external_seq_region_id, internal_seq_region_id FROM seq_region_mapping WHERE mapping_set_id=$max_mapping_set_id;");
  $sth->execute();
  while (my @row = $sth->fetchrow_array) {
    my $external_seq_region_id = $row[0]; # previous release
    my $internal_seq_region_id = $row[1]; # current release
    $id_mapping->{$external_seq_region_id} = $internal_seq_region_id;     
  }
  $sth->finish();

  my $vdbh = $variation_dba->dbc->db_handle;
  foreach my $prev_seq_region_id ( keys %$id_mapping) {
    my $new_seq_region_id = $id_mapping->{$prev_seq_region_id};
    if ($dry_run) {
      print "For $dbname: Update seq_region SET seq_region_id=$new_seq_region_id WHERE seq_region_id=$prev_seq_region_id\n";
    } else {
      $vdbh->do("Update seq_region SET seq_region_id=$new_seq_region_id WHERE seq_region_id=$prev_seq_region_id") or die $dbh->errstr;
    }
  }

}


sub get_max_mapping_set_id {
  my ($dbh, $dbname) = @_;
  my $sth_mapping = $dbh->prepare("select max(mapping_set_id) from $dbname.mapping_set");
  $sth_mapping->execute();
  my ($max_mapping_set_id) = $sth_mapping->fetchrow_array();
  if (!defined $max_mapping_set_id) { return 0; }
  return $max_mapping_set_id;
}


1;
