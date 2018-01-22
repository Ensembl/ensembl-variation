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


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk.org>.

=cut

use warnings;
use Getopt::Long;
use dbSNP::GenericChromosome;
use dbSNP::DBManager;

use vars qw(@ISA);
use ImportUtils qw(debug load dumpSQL create_and_load);
use Progress;

our($species,$shared_db,$registry_file);

GetOptions('species=s'         => \$species,
           'shared_db=s'       => \$shared_db,
           'registry_file=s'   => \$registry_file
          );

$species ||= 'human';

usage('-shared_db argument is required')     if (!$shared_db);
usage('-registry_file argument is required') if (!$registry_file);

my $db = "pontus_dbsnp132_human_external_data";

my $dbm = dbSNP::DBManager->new($registry_file,$species);
$dbm->dbSNP_shared($shared_db);
my $dbSNP  = $dbm->dbSNP()->dbc();
my $dbCore = $dbm->dbCore()->dbc();
my $dbVar  = $dbm->dbVar()->dbc();
  
my ($cs) = @{$dbm->dbCore()->get_CoordSystemAdaptor->fetch_all};
my $ASSEMBLY_VERSION = $cs->version();

# Main #
my $start = time();
dump_AFFYIDs();
print Progress::location();
my $end = time();
my $duration = Progress::time_format($end-$start);
print $duration->{'weeks'} . " weeks, " . $duration->{'days'} . " days, " . $duration->{'hours'} . " hours, " . $duration->{'minutes'} . " minutes and " . $duration->{'seconds'} . " seconds spent in dump_AFFYIDs()\n";

  
# This method uses the pontus_affy_array_mapping database on ens-variation
sub dump_AFFYIDs{

  my $self = shift;
  my ($source_name,$source_description,$source_url,$set_name);
    
  my $stmt;
  $source_url = "http://www.affymetrix.com/";
  
  debug("Dumping AFFY information from dbSNP");
  $stmt = qq{
    SELECT
  'rs'+LTRIM(STR(s.snp_id)),
  s.loc_snp_id,
  loc_batch_id
    FROM
  SubSNP s,
  Batch b
    WHERE
  s.batch_id = b.batch_id AND
  b.handle = 'AFFY'
  };
  dumpSQL($dbSNP,$stmt);
  
  debug("Loading  ids into temporary table");
  create_and_load($dbVar,"tmp_rs_AFFY","rsID *","AFFYid", "affy_name");
  print Progress::location();
  
  foreach my $table ("affy_array_name_pair_100k","affy_array_name_pair_500k","affy_array_name_pair_g6") {
    
    if ($table =~ /100k/i) {
      $source_name = "Affy GeneChip 100K Array";
      $source_description = "Variants from the Affymetrix GeneChip Human Mapping 100K Array Set";
      $set_name = "Mapping50K";
    }
    elsif ($table =~ /500k/i) {
      $source_name = "Affy GeneChip 500K Array";
      $source_description = "Variants from the Affymetrix GeneChip Human Mapping 500K Array Set";
      $set_name = "Mapping250K";
    }
    elsif ($table =~ /g6/i) {
      $source_name = "Affy GenomeWideSNP_6.0";
      $source_description = "Variants from the Affymetrix Genome-Wide Human SNP Array 6.0";
      $set_name = "6.0";
    }

    debug(localtime() . "\tCreating name_pair table with source_name $source_name...");
    $stmt = qq{
      CREATE TABLE $table LIKE $db\.$table
    };
    $dbVar->do($stmt);
    print Progress::location();
    
    debug(localtime() . "\tDumping AFFY id name pairs from $db\.$table\n");
    $stmt = qq{
      INSERT INTO $table ( affy_name, rs_name )
      SELECT affy_name, rs_name
      FROM $db\.$table
    };
    $dbVar->do($stmt);
    
    debug(localtime() . "\tCreating a temporary table for AFFY ids\n");
    $stmt = qq{ 
			 CREATE TABLE tmp
         SELECT
           a.AFFYid AS affy_name,
           a.rsID AS rs_name
         FROM
           tmp_rs_AFFY a,
           $table c
         WHERE
           a.AFFYid = c.affy_name
    };
    $dbVar->do($stmt);
    print Progress::location();
    
    debug(localtime() . "\tAdding AFFY ids from temporary table\n");
    $stmt = qq{
      INSERT IGNORE INTO
        $table (
          affy_name,
          rs_name
        )
        SELECT
          affy_name,
          rs_name
        FROM
          tmp
    };
    $dbVar->do($stmt);
    print Progress::location();
    
    debug(localtime() . "\tDropping temporary table\n");
    $stmt = qq{ DROP TABLE tmp };
    $dbVar->do($stmt);
    print Progress::location();
    
    
    # The code below has been replaced by the dump-and-load code above
    #$dbVar->do(qq{CREATE TABLE $table\_name_pair like $table.name_pair});
    #$dbVar->do(qq{insert into $table\_name_pair select * from $table.name_pair});
    #$dbVar->do(qq{insert ignore into $table\_name_pair
    #                        select a.AFFYid as affy_name,a.rsID as rs_name 
    #                        from tmp_rs_AFFY a, $table.name_pair c 
    #                        where a.AFFYid=c.affy_name});

    if ($table =~ /g6/i) {
      debug(localtime() . "\tUpdating $table\n");
      $stmt = qq{
        UPDATE $table SET affy_name = REPLACE(affy_name,'AFFY_6_1M_','')
      };
      $dbVar->do($stmt);
      print Progress::location();
    }
    
    my $source_id_ref = $dbVar->db_handle->selectall_arrayref(qq{
                       SELECT source_id from source where name = "$source_name"});
    my $source_id = $source_id_ref->[0][0];

    if (!$source_id) {
      $dbVar->do(qq{insert into source (name,description,url,type) 
		                values("$source_name","$source_description","$source_url","chip")
								   }
							  );
      print Progress::location();
      $source_id = $dbVar->db_handle->{'mysql_insertid'};
    }

    debug("Inserting in variation_synonym table from $table\...");
    $dbVar->do(qq{ INSERT IGNORE INTO variation_synonym (variation_id, source_id, name)
                      SELECT v.variation_id, $source_id as source_id, t.affy_name as name 
                      FROM variation v, $table t
                      WHERE v.name = t.rs_name
                 }
    );
    print Progress::location();

    #update rs_ID to rsCurrent from rsHigh
    #$dbVar->do(qq{update tmp_rs_AFFY_test t, rsHist h set t.rsID=h.rsCurrent where t.rsID=h.rsHigh});
    debug("Inserting in variation_synonym table from table  tmp_rs_AFFY with $source_name...");
    $dbVar->do(qq{ INSERT IGNORE INTO variation_synonym (variation_id, source_id, name)
                     SELECT v.variation_id, $source_id as source_id, t.AFFYid as name 
                     FROM variation v, tmp_rs_AFFY t
                     WHERE v.name = t.rsID
                       AND t.affy_name like "%$set_name%"
                 }
              );
    print Progress::location();

   #and finally, remove the temporary table
    $dbVar->do(qq{DROP TABLE $table});
    print Progress::location();

  }
}


sub usage {
    my $msg = shift;
    
    print STDERR <<EOF;
    
  usage: perl import_affy.pl <options>
      
  options:
    -species <species>              name of the species (default = human)
    -shared_db <shared db>          name of the shared DB (e.g. shared_134)
    -registry_file <registry file>  registry file name (e.g. './ensembl.registry')
      
EOF

      die("\n$msg\n\n");
}
