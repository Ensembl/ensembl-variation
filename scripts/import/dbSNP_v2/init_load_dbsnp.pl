#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2022] EMBL-European Bioinformatics Institute
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

## Script to set up seq_region table and tables for synonym_lookups

use strict;
use warnings;

use Getopt::Long;
use Bio::EnsEMBL::Registry;
use ImportUtils qw(debug load dumpSQL create_and_load);

my ($species, $registry_file, $tmpdir, $dbSNP_version);
my ($debug, $help);

GetOptions ("species=s"        => \$species,
            "registry=s"       => \$registry_file,
            "tmpdir=s"         => \$tmpdir,
            "dbsnp_version=i"  => \$dbSNP_version,
            "debug"            => \$debug,
            "help|h"           => \$help,
);

die usage() if ($help);

# Checking required parameters
die("A registry file (--registry) is required") unless (defined($registry_file));
die("A temporary directory (--tmpdir) is required") unless (defined($tmpdir));
die("The dbSNP version (--dbsnp_version) is required") unless (defined($dbSNP_version));

# Check the parameters are valid
if ($dbSNP_version !~ /^\d+$/) {
  die("The dbSNP_version ($dbSNP_version) is not numeric");
}

if (! -e $tmpdir) {
  die("The temporary directory ($tmpdir) does not exist");
}

if (! -e $registry_file) {
  die("The registry file ($registry_file) does not exist");
}

$species ||= 'homo_sapiens';

print "species = $species\n";
print "registry_file = $registry_file\n";
print "dbSNP_version = $dbSNP_version\n";
print "tmpdir = $tmpdir\n";

# Set up area to dump and load files from
$ImportUtils::TMP_DIR = $tmpdir;
$ImportUtils::TMP_FILE = 'dbtempfile';

my $reg = 'Bio::EnsEMBL::Registry';
$reg->no_version_check(1); 
$reg->load_all($registry_file);

my $dbVar = $reg->get_DBAdaptor($species, 'variation') || die "Error getting db adaptor variation\n";
my $dbh_var = $dbVar->dbc->db_handle;
get_db_info($dbh_var) if $debug;

my $dbCore = $reg->get_DBAdaptor($species, 'core') || die "Error getting db adaptor core\n";
my $dbh_core = $dbCore->dbc->db_handle;
get_db_info($dbh_core) if $debug;

insert_source($dbh_var, $dbSNP_version);

create_seq_region_table($dbh_core, $dbh_var);
create_coord_system_table($dbh_core, $dbh_var);

# TODO have single sub for nt/nw with parameter
#      remove code duplication
create_nt_synonyms($dbh_core, $dbh_var);
create_nw_synonyms($dbh_core, $dbh_var);
create_nc_synonyms($dbh_core, $dbh_var);

add_failed_variation_set($dbh_var);

sub insert_source {
  my ($dbh_var, $version) = @_;

  my $dbname = 'dbSNP';
  my $url = 'http://www.ncbi.nlm.nih.gov/projects/SNP/';
  my $somatic_status = 'mixed';
  
  my ($source_id, $name, $description, $data_types);

  my $sql = qq{INSERT INTO source
               (source_id, name, version, description, 
                 url, somatic_status, data_types)
               VALUES 
               (?, ?, $version, ?,
                "$url", "$somatic_status", ?)};
  
  my $sth = $dbh_var->prepare($sql);

  # Insert the dbSNP source
  $source_id = 1;
  $name = $dbname;
  $description = 'Variants (including SNPs and indels) imported from dbSNP';
  $data_types = 'variation';
  $sth->execute($source_id, $name, $description, $data_types);

  # Insert the Archive dbSNP - dbsnp v1 merges
  $source_id = 2;
  $name = "Archive $dbname";
  $description = 'Former dbSNP variant names, merged by variant';
  $data_types = 'variation_synonym';
  $sth->execute($source_id, $name, $description, $data_types);
  
  # Insert the dbSNP HGVS
  $source_id = 3;
  $name = "$dbname HGVS";
  $description = 'HGVS annotation from dbSNP';
  $sth->execute($source_id, $name, $description, $data_types);

  # Insert the Former dbSNP - dbsnp v2 merges 
  $source_id = 4;
  $name = "Former $dbname";
  $description = 'Former dbSNP variant names, merged by allele';
  $sth->execute($source_id, $name, $description, $data_types);
}

sub create_nw_synonyms {
  my ($dbh_core, $dbh_var) = @_;
  debug("\n>>>> create_nw_synonyms <<<<") if ($debug);
  if ($debug) {
    get_db_info($dbh_core);
    get_db_info($dbh_var);
  }
  
  dumpSQL($dbh_core, qq{
   SELECT sr2.seq_region_id, sr2.name, srs.synonym, assembly.asm_start, assembly.asm_end,assembly.ori
   FROM seq_region sr1, seq_region_synonym srs, seq_region sr2, assembly, seq_region_attrib sra, attrib_type at
   WHERE sr1.seq_region_id = srs.seq_region_id
   AND sr2.seq_region_id = assembly.asm_seq_region_id
   AND  sr1.seq_region_id = assembly.cmp_seq_region_id and srs.synonym like 'NW%'
   AND sra.attrib_type_id=at.attrib_type_id 
   AND at.code="toplevel" 
   AND sr2.seq_region_id = sra.seq_region_id}, 'MySQL');
   create_and_load($dbh_var, "tmp_nw_synonym", "seq_region_id i*", "name",
                                               "srs_synonym *", "asm_start i", 
                                               "asm_end i", "ori i");
}

sub create_nt_synonyms {
  my ($dbh_core, $dbh_var) = @_;
  debug("\n>>>> create_nt_synonyms <<<<") if ($debug);
  if ($debug) {
    get_db_info($dbh_core);
    get_db_info($dbh_var);
  }
  
  dumpSQL($dbh_core, qq{
   SELECT sr2.seq_region_id, sr2.name, srs.synonym, assembly.asm_start, assembly.asm_end,assembly.ori
   FROM seq_region sr1, seq_region_synonym srs, seq_region sr2, assembly, seq_region_attrib sra, attrib_type at
   WHERE sr1.seq_region_id = srs.seq_region_id
   AND sr2.seq_region_id = assembly.asm_seq_region_id
   AND  sr1.seq_region_id = assembly.cmp_seq_region_id and srs.synonym like 'NT%'
   AND sra.attrib_type_id=at.attrib_type_id 
   AND at.code="toplevel" 
   AND sr2.seq_region_id = sra.seq_region_id}, 'MySQL');
   create_and_load($dbh_var, "tmp_nt_synonym", "seq_region_id i*", "name",
                                               "srs_synonym *", "asm_start i", 
                                               "asm_end i", "ori i");
}

sub create_nc_synonyms {
  my ($dbh_core, $dbh_var) = @_;
  if ($debug) {
    debug("\n>>>> create_nc_synonyms <<<<") if ($debug);
    get_db_info($dbh_core);
    get_db_info($dbh_var);
  }

  dumpSQL($dbh_core, qq{
   SELECT s.seq_region_id, s.name, srs.synonym
   FROM seq_region_synonym srs, seq_region s, seq_region_attrib sra, attrib_type att
   WHERE srs.seq_region_id = s.seq_region_id
   AND s.seq_region_id = sra.seq_region_id
   AND sra.attrib_type_id = att.attrib_type_id
   AND att.code="toplevel"
   AND srs.synonym like 'NC%'}, 'MySQL');

   create_and_load($dbh_var, "tmp_nc_synonym", "seq_region_id i*", "name",
                                               "srs_synonym *");
}

sub create_seq_region_table {
  my ($dbh_core, $dbh_var) = @_;
  if ($debug) {
    debug("\n>>>> create_seq_region_table <<<<") if ($debug);
    get_db_info($dbh_core);
    get_db_info($dbh_var);
  }
  dumpSQL($dbh_core, qq{SELECT sr.seq_region_id, sr.name, sr.coord_system_id
                        FROM seq_region_attrib sra, attrib_type at, seq_region sr
                        WHERE sra.attrib_type_id=at.attrib_type_id 
                        AND at.code="toplevel" 
                        AND sr.seq_region_id = sra.seq_region_id 
                        }, 'MySQL');
  load($dbh_var, "seq_region", "seq_region_id", "name", "coord_system_id");
}

sub create_coord_system_table {
  my ($dbh_core, $dbh_var) = @_;
  if ($debug) {
    debug("\n>>>> create_coord_system_table <<<<") if ($debug);
    get_db_info($dbh_core);
    get_db_info($dbh_var);
  }
  dumpSQL($dbh_core, qq{SELECT coord_system_id, species_id, name, version, rank, attrib
                        FROM coord_system
                        }, 'MySQL');
  load($dbh_var, "coord_system", "coord_system_id", "species_id", "name", "version", "rank", "attrib");
}

sub add_failed_variation_set {
  my ($dbh_var) = @_;
  debug("\n>>>> add_failed_variation_set <<<<") if ($debug);

  my $fail_attrib_ext_sth  = $dbh_var->prepare(qq[ select at.attrib_id
                                                from attrib at, attrib_type att
                                                where att.code = 'short_name' 
                                                and att.attrib_type_id = at.attrib_type_id
                                                and at.value = 'fail_all'
                                            ]);
 
  my $variation_set_ext_sth  = $dbh_var->prepare(qq[ select variation_set_id
                                                        from variation_set
                                                         where name = ?
                                                         ]);

  my $variation_set_ins_sth  = $dbh_var->prepare(qq[ insert into variation_set
                                                        ( name, description, short_name_attrib_id)
                                                         values (?,?,?)
                                                      ]);
  ## check if already present
  $variation_set_ext_sth->execute('All failed variations')  || die "Failed to extract failed variant set id\n";
  my $failed_set_id = $variation_set_ext_sth->fetchall_arrayref();

  unless(defined $failed_set_id->[0]->[0]){
    ## no set entered - look up attrib for short name and enter set

    $fail_attrib_ext_sth->execute() || die "Failed to extract failed set attrib reasons\n";
    my $attrib = $fail_attrib_ext_sth->fetchall_arrayref();

    die "Exiting: Error - attribs not found. Load attribs then re-run\n" unless defined $attrib->[0]->[0] ;

    ## if attribs loaded, enter set
    $variation_set_ins_sth->execute( 'All failed variations',
                                     'Variations that have failed the Ensembl QC checks' ,
                                     $attrib->[0]->[0] )|| die "Failed to insert failed set\n"; 

    ## and pull out id to return
    $variation_set_ext_sth->execute('All failed variations')  || die "Failed to extract failed variant set id\n";
    $failed_set_id = $variation_set_ext_sth->fetchall_arrayref();      
  }

  return $failed_set_id->[0]->[0];
}

sub usage {

  print qq{
  Usage: perl init_load_dbsnp.pl [OPTION]

  Creates seq region lookup tables for dbsnp-v2 import.
  
  Populates the seq_region table and creates tables for synonym lookups (NC, NW, NT)

  Options:
  -h | --help          Displays this message and quit
  --species            Species defaults to homo_sapiens
  --registry           Registry file [required]
  --tmpdir             Temporary directory to write dump files [required]
  --dbsnp_version      dbSNP version being loaded (integer) [required]
  --debug              Flag for debugging
  }. "\n";
  exit(0);
}

sub get_db_info {
  my ($dbh) = @_;
  
  debug("\n>>>> get_db_info <<<<") if ($debug);
  
  my $sql = "SELECT now(), database(), user()";
  my $sth = $dbh->prepare($sql);
  $sth->execute();
  $sth->dump_results();
}
