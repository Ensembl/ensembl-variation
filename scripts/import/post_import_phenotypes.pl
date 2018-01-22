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

use strict;
use warnings;

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;

use DBI qw(:sql_types);
use FileHandle;
use Getopt::Long;

usage() if (!scalar(@ARGV));

my $config = {};

GetOptions(
  $config,
  'registry=s',
  'species=s',
  'clean_up',  
  'help!',
) or die "Error: Failed to parse command line arguments\n";

usage() if ($config->{help});

die ('A registry file is required (--registry)') unless (defined($config->{registry}));
die ('A species must be defiened (--species)') unless (defined($config->{species}));

main();

sub main {

  my $phenotype_attrib_id = '418';

  my $registry = 'Bio::EnsEMBL::Registry';   
  $registry->load_all($config->{registry});
  my $species = $config->{species};
  my $vdba = $registry->get_DBAdaptor($species, 'variation');
  my $dbh = $vdba->dbc->db_handle;

  if ($config->{clean_up}) {
    # clean up all evidence_attribs with previous phenotype attrib ids
    $dbh->do(qq{ DROP TABLE IF EXISTS variation_ids_old_phenotype_evdn});
    $dbh->do(qq{ CREATE TABLE `variation_ids_old_phenotype_evdn` (`variation_id` int(10) unsigned NOT NULL, PRIMARY KEY (`variation_id`))});
    $dbh->do(qq{ INSERT INTO variation_ids_old_phenotype_evdn SELECT distinct variation_id FROM variation WHERE evidence_attribs LIKE '%$phenotype_attrib_id%'}) ;
    $dbh->do(qq{
      UPDATE variation v, variation_ids_old_phenotype_evdn e
      SET v.evidence_attribs = REPLACE(v.evidence_attribs, '$phenotype_attrib_id', '')
      WHERE e.variation_id = v.variation_id;
    });
    $dbh->do(qq{
      UPDATE variation_feature vf, variation_ids_old_phenotype_evdn e
      SET vf.evidence_attribs = REPLACE(vf.evidence_attribs, '$phenotype_attrib_id', '')
      WHERE e.variation_id = vf.variation_id;
    });
  }
  # collect all variations with phenotype associations and update evidence_attribs
  $dbh->do(qq{ DROP TABLE IF EXISTS variation_ids_new_phenotype_evdn});
  $dbh->do(qq{ CREATE TABLE `variation_ids_new_phenotype_evdn` (`variation_id` int(10) unsigned NOT NULL, PRIMARY KEY (`variation_id`))});
  $dbh->do(qq{ INSERT INTO variation_ids_new_phenotype_evdn SELECT distinct v.variation_id FROM variation v, phenotype_feature pf WHERE pf.type = 'Variation' AND pf.object_id = v.name; }) ;

  $dbh->do(qq{
    UPDATE variation v, variation_ids_new_phenotype_evdn e
    SET v.evidence_attribs = CONCAT_WS(',', v.evidence_attribs, '$phenotype_attrib_id')
    WHERE e.variation_id = v.variation_id;
  });
  $dbh->do(qq{
    UPDATE variation_feature vf, variation_ids_new_phenotype_evdn e
    SET vf.evidence_attribs = CONCAT_WS(',', vf.evidence_attribs, '$phenotype_attrib_id')
    WHERE e.variation_id = vf.variation_id;
  });

  # update display column
  $dbh->do(qq{
    UPDATE variation v, variation_ids_new_phenotype_evdn e
    SET v.display = 1
    WHERE e.variation_id = v.variation_id;
  });
  $dbh->do(qq{
    UPDATE variation_feature vf, variation_ids_new_phenotype_evdn e
    SET vf.display = 1
    WHERE e.variation_id = vf.variation_id;
  });
  $dbh->do(qq{
    UPDATE transcript_variation tv, variation_feature vf, variation_ids_new_phenotype_evdn e
    SET tv.display = 1
    WHERE e.variation_id = vf.variation_id
    AND vf.variation_feature_id = tv.variation_feature_id;
  });
  $dbh->do(qq{
    UPDATE MTMP_transcript_variation tv, variation_feature vf, variation_ids_new_phenotype_evdn e
    SET tv.display = 1
    WHERE e.variation_id = vf.variation_id
    AND vf.variation_feature_id = tv.variation_feature_id;
  });


}

sub usage {

  print qq{
  Usage: perl post_import_phenotype_data.pl -registry [registry_file] -species [species] [OPTIONS]

  Update evidence_attribs column in variation and variation_feature tables for all variants with
  phenotype asssociations.
  Set display to 1 in variation and variation_feature if variation has phenotype associations.

  Options:
    -help        Print this message
    -clean_up    Use to clean up all Phenotype_or_Disease evidence attributes from the previous release
  } . "\n";
  exit(0);
}
