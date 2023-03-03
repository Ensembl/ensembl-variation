=head1 LICENSE
Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2023] EMBL-European Bioinformatics Institute
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
use Bio::EnsEMBL::Variation::Utils::Date;
use POSIX;

sub fetch_input {
  my $self = shift;
}

sub run {
  my $self = shift;

  my %tables = (
    "batch_variation" => ["name", "source_id",
                          "minor_allele", "minor_allele_freq", "minor_allele_count",
                          "evidence_attribs", "display", "class_attrib_id"],
    "batch" => ["filename", "parent_filename"],
    "failed_variation_feature_spdi" => ["variation_feature_id", "spdi_failed_description_id"],
    "failed_variation_feature" => ["variation_feature_id", "failed_description_id"],
    "failed_variation" => ["variation_id", "failed_description_id"],
    "placement_allele" => ["variation_id", "seq_id", "position",
                           "deleted_sequence", "inserted_sequence", "hgvs"],
    "tmp_variation_citation" => ["variation_id", "pmid"],
    "variation" => ["variation_id", "name", "source_id",
                    "evidence_attribs", "display", "class_attrib_id"],
    "variation_feature" =>["variation_feature_id" ,"variation_name", "map_weight", "seq_region_id",
                           "seq_region_start", "seq_region_end", "seq_region_strand",
                           "variation_id", "allele_string", "ancestral_allele",
                           "source_id", "variation_set_id", "class_attrib_id",
                           "evidence_attribs", "display"],
    "variation_synonym" => ["variation_id", "source_id", "name"],
  );

  # The input_directory is <data-dir>/<sub_dir>

  $ImportUtils::TMP_DIR = $TMP_DIR;
  $ImportUtils::TMP_FILE = "variation_" . $TMP_FILE;
  ## Load variation dump file
  load($config->{'dbh_var'},
    ("variation", 
      ("variation_id", "name", "source_id", "evidence_attribs",
      "display", "class_attrib_id")
    )
  );

  $ImportUtils::TMP_FILE = "variation_feature_" . $TMP_FILE;
  # ## Load variation_feature dump file
  load($config->{'dbh_var'},
    ("variation_feature", 
      ("variation_name", "map_weight", "seq_region_id",
      "seq_region_start", "seq_region_end", "seq_region_strand",
      "variation_id", "allele_string", "ancestral_allele",
      "source_id", "variation_set_id", "class_attrib_id",
      "evidence_attribs", "display")
    )
  );  

}

1;
