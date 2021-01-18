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

package Bio::EnsEMBL::Variation::Pipeline::Remapping::ComparePreviousAssembly;

use base ('Bio::EnsEMBL::Variation::Pipeline::Remapping::BaseRemapping');

use Bio::EnsEMBL::Registry;
use FileHandle;
use List::MoreUtils qw(uniq);
sub fetch_input {
  my $self = shift;
  $self->SUPER::fetch_input;
}

sub run {
  my $self = shift;
  my $feature_table = $self->param('feature_table');
  my $vdba_oldasm = $self->get_oldasm_variation_database_connection;
  my $counts_prev_assembly = $self->get_feature_counts_by_seq_region($vdba_oldasm, $feature_table);
  my $vdba_newasm = $self->get_newasm_variation_database_connection;
  my $counts_new_assembly = $self->get_feature_counts_by_seq_region($vdba_newasm, $feature_table);
  my $pipeline_dir = $self->param('pipeline_dir');
  my $fh = FileHandle->new("$pipeline_dir/QC_ComparePreviousAssembly", 'w');
  my @seq_regions = uniq (keys %$counts_prev_assembly, keys %$counts_new_assembly);

  my $count_old2new_loss = 0;
  my $count_old2new_gain = 0;
  my $count_suspicious_decrease = 0;
  my $count_suspicious_increase = 0;

  print $fh "seq_region count_prev count_new ratio\n";
  foreach my $seq_region (@seq_regions) {
    my $count_prev = $counts_prev_assembly->{$seq_region} || 'NA';
    my $count_new =  $counts_new_assembly->{$seq_region} || 'NA';
    $count_old2new_loss++ if ($count_new eq 'NA');
    $count_old2new_gain++ if ($count_prev eq 'NA');
    if ("$count_prev" ne 'NA' && "$count_new" ne 'NA') {
      my $ratio = $count_new / $count_prev;
      $count_suspicious_increase++ if ($ratio > 1.2);
      $count_suspicious_decrease++ if ($ratio < 0.8);
      print $fh "$seq_region $count_prev $count_new $ratio\n";
    } else {
      print $fh "$seq_region $count_prev $count_new\n";
    }
  } 
  $fh->close; 

  $fh = FileHandle->new("$pipeline_dir/QC_ComparePreviousAssembly_report", 'w');
  print $fh "Lost variants on $count_old2new_loss seq_regions\n" if ($count_old2new_loss);
  print $fh "Gained variants on $count_old2new_gain seq_regions\n" if ($count_old2new_gain);
  print $fh "$count_suspicious_decrease seq_regions with suspiciously high decrease in number of variants\n" if ($count_suspicious_decrease);
  print $fh "$count_suspicious_increase seq_regions with suspiciously high increase in number of variants\n" if ($count_suspicious_increase);

  $fh->close;
}

sub get_feature_counts_by_seq_region {
  my $self = shift;
  my $vdba = shift;
  my $feature_table = shift;
  my $dbh = $vdba->dbc->db_handle();
  my $extra_sql = '';
  if ($self->param('mode') eq 'remap_QTL') {
    $extra_sql = " WHERE ft.type='QTL' ";
  }
  my $sth = $dbh->prepare(qq{
    SELECT sr.name, count(ft.seq_region_id)
    FROM $feature_table ft
    LEFT JOIN
    seq_region sr
    ON (ft.seq_region_id = sr.seq_region_id)
    $extra_sql
    GROUP BY sr.name;
  }, {mysql_use_result => 1});

  $sth->execute() or die $sth->errstr;
  my $seq_region_2_count = {};
  while (my $row = $sth->fetchrow_arrayref) {
    my ($seq_region, $count) = @$row;
    $seq_region_2_count->{$seq_region} = $count;
  }
  $sth->finish();
  return $seq_region_2_count;
}

1;
