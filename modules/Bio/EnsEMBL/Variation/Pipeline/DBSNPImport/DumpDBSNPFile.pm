=head1 LICENSE
Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2026] EMBL-European Bioinformatics Institute
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

=head1 NAME

Bio::EnsEMBL::Variation::Pipeline::DBSNPImport::DumpDBSNPFile

=head1 DESCRIPTION

Dump a dbSNP tmp files into databases

=cut

package Bio::EnsEMBL::Variation::Pipeline::DBSNPImport::DumpDBSNPFile;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);

use FileHandle;
use Bio::EnsEMBL::Registry;
use File::Path qw(make_path);
use ImportUtils qw(load);
use List::MoreUtils qw(first_index);
use File::Basename;
use Bio::EnsEMBL::Variation::Utils::Date;
use POSIX;

my @chrs = (
  "chr1", "chr2", "chr3", "chr4", "chr5",
  "chr6", "chr7", "chr8", "chr9", "chr10",
  "chr11", "chr12", "chr13", "chr14", "chr15",
  "chr16", "chr17", "chr18", "chr19", "chr20",
  "chr21", "chr22", "chrX" , "chrY", "chrMT",
  "other"
);

my %ids = (
  "variation_id" => 0,
  "variation_feature_id"  => 0,
  "batch_id"  => 0,
  "failed_variation_id" => 0,
  "failed_variation_feature_id" => 0,
  "placement_allele_id"  => 0,
  "variation_synonym_id"  => 0
);

my %tables = (
  "batch_variation" => ["batch_id", "variation_id", "variant_type"],
  "batch" => ["batch_id", "filename", "parent_filename"],
  "failed_variation_feature_spdi" => ["variation_feature_id", "spdi_failed_description_id"],
  "failed_variation_feature" => ["failed_variation_feature_id", "variation_feature_id", "failed_description_id"],
  "failed_variation" => ["failed_variation_id", "variation_id", "failed_description_id"],
  "placement_allele" => ["placement_allele_id", "variation_id", "seq_id", "position",
                          "deleted_sequence", "inserted_sequence", "hgvs"],
  "tmp_variation_citation" => ["variation_id", "pmid"],
  "variation" => ["variation_id", "name", "source_id",
                  "evidence_attribs", "display", "class_attrib_id"],
  "variation_feature" =>["variation_feature_id" ,"variation_name", "map_weight", "seq_region_id",
                          "seq_region_start", "seq_region_end", "seq_region_strand",
                          "variation_id", "allele_string", "ancestral_allele",
                          "source_id", "variation_set_id", "class_attrib_id",
                          "evidence_attribs", "display"],
  "variation_synonym" => ["variation_synonym_id", "variation_id", "source_id", "name"],
);

sub fetch_input {
  my $self = shift;
}

sub run {
  my $self = shift;

  my $var_dba = $self->get_species_adaptor('variation');
  my $dbh = $var_dba->dbc();

  my $TMP_DIR = $self->param_required("pipeline_dir");

  # The input_directory is <data-dir>/<sub_dir>

  for my $chr (@chrs) {
    my @return_value = glob("ls ${TMP_DIR}/split-src/${chr}/tmp*");
    my @chr_files = grep( /${chr}/, @return_value);
    for my $file (@chr_files) {
      my $base = basename($file);
      $base =~ m/tmp_(.*)_refsnp/;
      my $table = $1;
      
      my $count = $ids{"${table}_id"};

      if (-e "${TMP_DIR}/tmp_${table}.txt"){
        $self->increment_and_load_file(${table}, $file, "${TMP_DIR}/tmp_${table}.txt");
        $ids{"${table}_id"} = count_lines("$file", $count) if (defined($count));
      } else {
        $self->run_system_command("cat ${file} > ${TMP_DIR}/tmp_${table}.txt");
        $ids{"${table}_id"} = count_lines("$file", $count) if (defined($count));
      }
    }
  }

  ## Sort variation_feature table
  my ($return_value, $stderr, $flat_cmd) = $self->run_system_command("sort --parallel=8 --buffer-size=6G -k4,4n -k5,5n -k6,6n -o ${TMP_DIR}/sorted_tmp_variation_feature.txt ${TMP_DIR}/tmp_variation_feature.txt");
  if ($return_value) {
      die("there was an error running as ($flat_cmd: $stderr)");
  } else {
    $self->run_system_command("rm ${TMP_DIR}/tmp_variation_feature.txt");
    $self->run_system_command("mv ${TMP_DIR}/sorted_tmp_variation_feature.txt ${TMP_DIR}/tmp_variation_feature.txt");
  }

  ## For loop to load all sorted tables
  for my $table (keys %tables) {
    $ImportUtils::TMP_DIR = $TMP_DIR;
    $ImportUtils::TMP_FILE = "tmp_${table}.txt";
    ## Load variation dump file
    load($dbh,
      ("${table}", 
        @{$tables{"${table}"}}
      )
    );
  }

}

sub count_lines {
  my ($file, $count) = @_;

  open(FILE, "< $file") or die "can't open $file: $!";
  $count += tr/\n/\n/ while sysread(FILE, $_, 2 ** 16);

  return $count;

}

sub increment_and_load_file {
  my $self = shift;
  my $table = shift;
  my $file = shift;
  my $tmp = shift;

  # Find all col_idx
  my ($cmd_string, $count);

  for my $id (keys %ids) {
    my $col_idx = first_index { $_ eq $id } @{$tables{"$table"}};
    next if $col_idx == -1;
    $col_idx++;

    $count = $ids{"$id"};

    my $sub_string = "\\\$${col_idx}=\\\$${col_idx}+${count}";
    $cmd_string .= defined($cmd_string)? ";${sub_string}" : "${sub_string}";
  }

  my $full_cmd = 'awk -F' . '"\\t" ' . '"{' . ${cmd_string} . ';print}" OFS="\\t" ' . ${file} . ' >> ' . ${tmp};

  # Increment + add to main tmp
  $self->run_system_command($full_cmd);

}

1;
