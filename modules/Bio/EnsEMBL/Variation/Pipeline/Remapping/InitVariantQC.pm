#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

package Bio::EnsEMBL::Variation::Pipeline::Remapping::InitVariantQC;

use strict;
use warnings;

use FileHandle;
use Bio::EnsEMBL::Registry;
use base ('Bio::EnsEMBL::Variation::Pipeline::Remapping::BaseRemapping');

sub fetch_input {
  my $self = shift;
  my $registry = 'Bio::EnsEMBL::Registry';
  $registry->load_all($self->param('registry_file_newasm'));  
  my $vdba = $registry->get_DBAdaptor($self->param('species'), 'variation');
  $self->param('vdba', $vdba);
}

sub run {
  my $self = shift;
  my $vdba = $self->param('vdba');

  my $dbh = $vdba->dbc->db_handle();
  my $pipeline_dir = $self->param('pipeline_dir');

  my $dump_mapped_features_dir = $self->param('dump_mapped_features_dir');
  
  my $file_count = 1;
  my $entries_per_file = $self->param('entries_per_file'); 
  my $count_entries = 0;
  my $fh = FileHandle->new("$dump_mapped_features_dir/$file_count.txt", 'w');

  my @column_names = qw/variation_feature_id variation_id seq_region_id seq_region_start seq_region_end seq_region_strand allele_string map_weight/;
  my $column_concat = join(',', @column_names);
 
  my $feature_table = $self->param('feature_table');
  if ($self->param('mode') eq 'remap_post_projection') {
    $feature_table = "$feature_table\_post_projection";
  }
 
  my $sth = $dbh->prepare(qq{SELECT $column_concat FROM $feature_table}, {mysql_use_result => 1});
        
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
      $fh = FileHandle->new("$dump_mapped_features_dir/$file_count.txt", 'w');
      $count_entries = 0;
    } 
    $count_entries++;
    print $fh join("\t", @pairs), "\n";
  }

  $sth->finish();
  $fh->close();
  $self->param('file_count', $file_count);

  # create backup of failed_variation table
  my $failed_variation_table_BU = 'failed_variation_pre_remapping';
  $dbh->do(qq{ DROP TABLE IF EXISTS $failed_variation_table_BU});
  $dbh->do(qq{ CREATE TABLE $failed_variation_table_BU like failed_variation;});
  $dbh->do(qq{ INSERT INTO $failed_variation_table_BU select * FROM failed_variation;});

  # update failed_variants map_weight > 1 
  # Variation maps to more than one genomic location
  my $failed_description_id = 19;
  $dbh->do(qq{DELETE FROM failed_variation WHERE failed_description_id=$failed_description_id;});
  $dbh->do(qq{INSERT INTO failed_variation(variation_id, failed_description_id) SELECT distinct variation_id, $failed_description_id from $feature_table WHERE map_weight > 1;});

  # update failed_variants no mapping  
  # Variation can not be re-mapped to the current assembly
  $failed_description_id = 17; 
  $dbh->do(qq{DELETE FROM failed_variation WHERE failed_description_id=$failed_description_id;});
  my $unmapped_variations = {};
  my $filtered_mapping_dir = $self->param('filtered_mappings_dir'); 
  opendir(DIR, $filtered_mapping_dir) or die $!;
  while (my $file = readdir(DIR)) {
    if ($file =~ /^failed/) {
      my $fh = FileHandle->new("$filtered_mapping_dir/$file", 'r');
      while (<$fh>) {
        chomp;
        my ($query_name, $type) = split /\t/;
        # query_name: 156358-150-1-150-11:5502587:5502587:1:T/C:rs202026261:dbSNP:SNV
        my @query_name_components = split('-', $query_name, 2);
        my $snd_part = $query_name_components[1];
        my $variation_name = (split(':', $snd_part))[5];
        $unmapped_variations->{"'$variation_name'"} = 1;
      }
      $fh->close();
    }
  }
  closedir(DIR);
  my $variation_names_concat = join(',', keys %$unmapped_variations);
  $dbh->do(qq{INSERT INTO failed_variation(variation_id, failed_description_id) SELECT variation_id, $failed_description_id FROM variation WHERE name IN ($variation_names_concat);});

  # update display column in variation and variation_feature tables

  # evidence_attribs: 371 = Cited
  $dbh->do(qq{INSERT INTO failed_variation(variation_id, failed_description_id) SELECT variation_id, $failed_description_id FROM variation WHERE name IN ($variation_names_concat);});

  # reset display column, set all values to 1
  $dbh->do(qq{UPDATE variation SET display=1 WHERE display=0;});
  $dbh->do(qq{UPDATE $feature_table SET display=1 WHERE display=0;});
  # set display to 0 for failed_variations
  $dbh->do(qq{UPDATE variation v, failed_variation fv SET v.display=0 WHERE fv.variation_id=v.variation_id;});
  $dbh->do(qq{UPDATE $feature_table rt, failed_variation fv SET rt.display=0 WHERE fv.variation_id=rt.variation_id;});
  # set display to 1 if variation is cited 
  $dbh->do(qq{UPDATE variation SET display=1 WHERE evidence_attribs LIKE '%371%';});
  $dbh->do(qq{UPDATE $feature_table rt, variation v SET rt.display=1 WHERE v.variation_id = rt.variation_id AND v.evidence_attribs LIKE '%371%';});
#ALTER TABLE variation_feature DROP COLUMN variation_feature_id_old 
#alter table variation_feature order by seq_region_id,seq_region_start,seq_region_end  
}

=begin
sub write_output {
  my $self = shift;
  my @jobs = ();
  my $file_count = $self->param('file_count');
  my $i = 1;
  while ($i <= $file_count) {
    push @jobs, {
      'file_number' => $i,
    };
    $i++;
  } 
  $self->dataflow_output_id(\@jobs, 2)
}
=end
=cut


1;
