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


use strict;
use warnings;

use Getopt::Long;
use ImportUtils qw(debug load dumpSQL);
use Bio::EnsEMBL::Utils::Exception qw(warning throw verbose);
use DBH;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use FindBin qw( $Bin );
use DBI qw(:sql_types);
use Data::Dumper;

my $species;

GetOptions('tmpdir=s'  => \$ImportUtils::TMP_DIR,
	   'tmpfile=s' => \$ImportUtils::TMP_FILE,
	   'species=s' => \$species
	   );

warn("Make sure you have a updated ensembl.registry file!\n");
die "you must specify the species !! " if (!$species);
my $registry_file ||= $Bin . "/ensembl.registry";

Bio::EnsEMBL::Registry->load_all( $registry_file );

my $dbVar = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'variation');
my $dbSanger = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'sanger');

my $TMP_DIR  = $ImportUtils::TMP_DIR;
my $TMP_FILE = $ImportUtils::TMP_FILE;

print "Time starting to import data: ", scalar(localtime),"\n";
my $call = "bsub -J import_job -o output_import_snps.txt -m bc_hosts/usr/bin/env perl import_Sanger_database.pl -tmpdir $TMP_DIR -tmpfile $TMP_FILE -species $species";
system($call);
sleep(10);
$call = "bsub -K -w 'done(import_job)' -J waiting_process sleep 1";
system($call);

rename "$TMP_DIR/variation.txt","$TMP_DIR/$TMP_FILE";
load($dbVar->dbc,"variation", "variation_id", "source_id", "name", "validation_status", "ancestral_allele");
 rename "$TMP_DIR/allele.txt","$TMP_DIR/$TMP_FILE";
 load($dbVar->dbc,"allele", "variation_id", "sample_id", "allele", "frequency");
 rename "$TMP_DIR/flanking_sequence.txt","$TMP_DIR/$TMP_FILE";
 load($dbVar->dbc,"flanking_sequence", "variation_id", "up_seq", "down_seq", "up_seq_region_start", "up_seq_region_end", "down_seq_region_start", "down_seq_region_end", "seq_region_id", "seq_region_strand");
# rename "$TMP_DIR/variation_synonym.txt","$TMP_DIR/$TMP_FILE";
# load($dbVar->dbc,"variation_synonym", "variation_id", "source_id", "name", "moltype");
 rename "$TMP_DIR/variation_feature.txt","$TMP_DIR/$TMP_FILE";
 load($dbVar->dbc,"variation_feature", "variation_feature_id", "variation_id", "seq_region_id", "seq_region_start", "seq_region_end", "seq_region_strand", "allele_string", "variation_name", "map_weight", "flags", "source_id", "validation_status", "consequence_type");
 rename "$TMP_DIR/transcript_variation.txt","$TMP_DIR/$TMP_FILE";
 load($dbVar->dbc,"transcript_variation", "transcript_id", "variation_feature_id", "cdna_start", "cdna_end", "translation_start", "translation_end", "peptide_allele_string", "consequence_type");
 rename "$TMP_DIR/read_coverage.txt","$TMP_DIR/$TMP_FILE";
 load($dbVar->dbc,"read_coverage", "seq_region_id", "seq_region_start", "seq_region_end", "level", "sample_id");
 rename "$TMP_DIR/tmp_individual_genotype_single_bp.txt","$TMP_DIR/$TMP_FILE";
 load($dbVar->dbc,"tmp_individual_genotype_single_bp","variation_id", "sample_id", "allele_1", "allele_2");

print "Time finished importing data: ", scalar(localtime),"\n";
