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
package Bio::EnsEMBL::Variation::Pipeline::RemappingVCF::InitRemapping;

use strict;
use warnings;

use FileHandle;
use Bio::EnsEMBL::Registry;

use base qw(Bio::EnsEMBL::Hive::Process);

sub run {
  my $self = shift;
  my $registry_file = $self->param('registry_file_newasm');
  my $species      = $self->param('species');
  my $vcf_file = $self->param('vcf_file');
  my $output = `tabix --list-chroms $vcf_file`;
  my @chroms_list = split/\n/, $output;


  my $registry = 'Bio::EnsEMBL::Registry';
  $registry->load_all($registry_file);
  my $dbh = $registry->get_DBAdaptor($species, 'variation')->dbc->db_handle;
  my @input = ();

  foreach my $vcf_chrom (@chroms_list) {
    my $chrom = $vcf_chrom;
    $chrom =~ s/chr//;
    if ($self->has_mappings($dbh, $chrom)) {
      push @input, {
        chrom => $chrom,
        vcf_chrom => $vcf_chrom,
      };
    }
  }

  $registry->clear;
  my $registry_file_oldasm = $self->param('registry_file_oldasm');
  $registry->load_all($registry_file_oldasm, 0, 1);
  my $slice_adaptor = $registry->get_adaptor($species, 'core', 'slice');

  my @trimmed_input = ();

  foreach my $input_hash (@input) {
    my $chrom = $input_hash->{chrom};
    my $slice = $slice_adaptor->fetch_by_region('toplevel', $chrom);
    if ($slice) {
      my $chrom_end = $slice->seq_region_end;
      $input_hash->{seq_region_end} = $chrom_end;
      push @trimmed_input, $input_hash;
    } else {
      $self->warning("No slice for $chrom");
    }
  }

  $self->param('input', \@trimmed_input);
}

sub write_output {
  my $self = shift;
  $self->dataflow_output_id($self->param('input'), 2);
}

sub has_mappings {
  my $self = shift;
  my $dbh = shift;
  my $seq_region_name = shift;

  my $sth = $dbh->prepare(qq{
    SELECT variation_id FROM vcf_variation
    WHERE seq_region_name_old = ?
    AND variation_id IS NOT NULL
    LIMIT 1;
  }, {mysql_use_result => 1});
  $sth->execute($seq_region_name);
  my @row = $sth->fetchrow_array;
  $sth->finish();
  return $row[0];
}

1;
