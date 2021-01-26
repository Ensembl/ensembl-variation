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
package Bio::EnsEMBL::Variation::Pipeline::RemappingVCF::LoadFromVCF;

use strict;
use warnings;

use FileHandle;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::Utils::Sequence qw(trim_sequences);
use Bio::EnsEMBL::IO::Parser::VCF4Tabix;
use Bio::EnsEMBL::IO::Parser::VCF4;
use ImportUtils qw(load);

use base qw(Bio::EnsEMBL::Hive::Process);

sub run {
  my $self = shift;
  $self->dump_data_from_VCF if ($self->param('dump_data_from_VCF'));
  $self->load_data if ($self->param('load_data'));
  $self->update_mappings if ($self->param('update_mappings'));
}

sub write_output {
  my $self = shift;
  $self->dataflow_output_id($self->param('input'), 2);
}

sub dump_data_from_VCF {
  my $self = shift;
  my $pipeline_dir = $self->param('pipeline_dir');
  my $registry_file = $self->param('registry_file_oldasm');
  my $species      = $self->param('species');
  my $registry = 'Bio::EnsEMBL::Registry';
  $registry->load_all($registry_file);
  my $slice_adaptor = $registry->get_adaptor($species, 'core', 'slice');
  my $vcf_file = $self->param('vcf_file');
  my $output = `tabix --list-chroms $vcf_file`;
  my @chroms = split/\n/, $output;

  my $fh = FileHandle->new("$pipeline_dir/vcf_variants_old_assembly", 'w');

  my $parser = Bio::EnsEMBL::IO::Parser::VCF4Tabix->open($vcf_file);

  foreach my $chrom (@chroms) {
    my $ensembl_chrom = $chrom;
    $ensembl_chrom =~ s/chr//;
    my $slice = $slice_adaptor->fetch_by_region('toplevel', $ensembl_chrom);
    my $seq_region_id_old = $slice->get_seq_region_id;
    my $chrom_end = $slice->seq_region_end;
    $parser->seek($chrom, 1, $chrom_end);
    while ($parser->next) {
      my $seq_name_old = $parser->get_seqname;
      my $start_old = $parser->get_start;
      my $start_padded_old = $parser->get_raw_start;
      my $reference = $parser->get_reference;
      my @alternatives = @{$parser->get_alternatives};
      my $allele_string = join('/', $reference, @alternatives);
      my $allele_string_padded = '';
      if ($start_old != $start_padded_old) {
        $allele_string_padded = $allele_string;
        my @trimmed_alleles = ();
        my $trimmed_ref = '';
        foreach my $allele (@alternatives) {
          my ($_ref, $_alt, $_start) = @{trim_sequences($reference, $allele, 1, 1, 1)};
          $trimmed_ref = $_ref;
          push @trimmed_alleles, $_alt;
        }
        $allele_string = join('/', $trimmed_ref, @trimmed_alleles);
      }
      print $fh join("\t", $ensembl_chrom, $seq_region_id_old, $start_old, $start_padded_old, $allele_string, $allele_string_padded), "\n";
    }
  }
  $fh->close;
}


sub load_data {
  my $self = shift; 
  my $pipeline_dir = $self->param('pipeline_dir');
  my $registry_file = $self->param('registry_file_newasm');
  my $species = $self->param('species');
  my $registry = 'Bio::EnsEMBL::Registry';
  $registry->load_all($registry_file);

  my $dbc = $registry->get_DBAdaptor($species, 'variation')->dbc;

  $dbc->do(qq{ DROP TABLE IF EXISTS vcf_variation}) or die $!;

  $dbc->do(qq{
    CREATE TABLE vcf_variation (
      vcf_variation_id int(10) unsigned NOT NULL AUTO_INCREMENT,
      seq_region_name_old varchar(50) DEFAULT NULL,
      seq_region_id_old int(10) unsigned DEFAULT NULL,
      seq_region_start_old int(11) DEFAULT NULL,
      seq_region_start_padded_old int(11) DEFAULT NULL,
      allele_string_old varchar(50) DEFAULT NULL,
      allele_string_padded_old varchar(50) DEFAULT NULL,
      vcf_id int(10) unsigned DEFAULT NULL,
      variation_id int(10) unsigned DEFAULT NULL,
      seq_region_name_new varchar(50) DEFAULT NULL,
      seq_region_id_new int(10) unsigned DEFAULT NULL,
      seq_region_start_new int(11) DEFAULT NULL,
      seq_region_start_padded_new int(11) DEFAULT NULL,
      allele_string_new varchar(50) DEFAULT NULL,
      allele_string_padded_new varchar(50) DEFAULT NULL,
      PRIMARY KEY (vcf_variation_id),
      KEY variation_idx (variation_id),
      KEY seq_region_name_old_idx (seq_region_name_old),
      KEY seq_region_id_old_idx (seq_region_id_old),
      KEY seq_region_start_old_idx (seq_region_start_old)
    );
  }) or die $!;


  my $dbh = $dbc->db_handle;
  
  my $TMP_DIR = $pipeline_dir;
  my $tmp_file = 'vcf_variants_old_assembly';
  my $result_table = 'vcf_variation';
  my $column_names_concat = 'seq_region_name_old,seq_region_id_old, seq_region_start_old,seq_region_start_padded_old,allele_string_old,allele_string_padded_old';
  $ImportUtils::TMP_DIR = $TMP_DIR;
  $ImportUtils::TMP_FILE = $tmp_file;
  load($dbh, ($result_table, $column_names_concat));
}

sub update_mappings {
  my $self = shift;

  my $registry = 'Bio::EnsEMBL::Registry';
  my $registry_file = $self->param('registry_file_oldasm_same_server');
  $registry->load_all($registry_file);
  my $species      = $self->param('species');
  my $dbc = $registry->get_DBAdaptor($species, 'variation')->dbc;
  my $dbname_oldasm = $dbc->dbname;

  $registry->clear;
  my $registry_file_newasm = $self->param('registry_file_newasm');
  $registry->load_all($registry_file_newasm, 0, 1);
  my $dbc = $registry->get_DBAdaptor($species, 'variation')->dbc;

  my $rows = $dbc->do(qq{
   UPDATE vcf_variation vcf, $dbname_oldasm.variation_feature vf
   SET vcf.variation_id = vf.variation_id
   WHERE vcf.seq_region_id_old = vf.seq_region_id
   AND vcf.seq_region_start_old = vf.seq_region_start;
  }) or die $dbc->errstr;

  my $vcf_variation_count = $self->_run_sql_query(qq{SELECT count(vcf_variation_id) FROM vcf_variation}, $dbc);
  my $variation_id_count = $self->_run_sql_query(qq{SELECT count(variation_id) FROM vcf_variation WHERE variation_id IS NOT NULL}, $dbc);
  my $ratio = $variation_id_count / $vcf_variation_count; 
  die "suspiciously low number of variation ids: $ratio. Rows affected $rows" if ($ratio < 0.8);

  $dbc->do(qq{
    UPDATE vcf_variation vcf, variation_feature vf, seq_region sr
    SET vcf.seq_region_name_new = sr.name, vcf.seq_region_id_new = vf.seq_region_id, vcf.seq_region_start_new = vf.seq_region_start, vcf.allele_string_new = vf.allele_string
    WHERE vcf.variation_id = vf.variation_id
    AND vf.seq_region_id = sr.seq_region_id;
  }) or die $dbc->errstr;

  my $seq_region_id_new_count = $self->_run_sql_query(qq{SELECT count(seq_region_id_new) FROM vcf_variation WHERE seq_region_id_new IS NOT NULL}, $dbc);
  my $ratio = $seq_region_id_new_count / $vcf_variation_count; 
  die "suspiciously low number of new seq region ids: $ratio" if ($ratio < 0.8);

}

sub _run_sql_query {
  my $self = shift;
  my $query = shift;
  my $dbc = shift;
  my $dbh = $dbc->db_handle;
  my $sth = $dbh->prepare($query);
  $sth->execute();
  my $row = $sth->fetchrow_arrayref;
  $sth->finish;
  return $row->[0];
}

sub write_output {
  my $self = shift;
}

1;
