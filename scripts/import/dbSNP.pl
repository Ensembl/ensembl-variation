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

# this is the script to fill the new variation
# schema with data from dbSNP
# we use the local mssql copy of dbSNP at the sanger center
# the script will call the dbSNP factory to create the object that will deal with
# the creation of the data according to the species


use strict;
use warnings;
use DBI;
use Getopt::Long;
use Bio::EnsEMBL::Registry;
use FindBin qw( $Bin );
use Progress;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw);

use ImportUtils qw(dumpSQL debug create_and_load load);
use dbSNP::GenericContig;
use dbSNP::GenericChromosome;
use dbSNP::MappingChromosome;
use dbSNP::Mosquito;
use dbSNP::Human;
use dbSNP::EnsemblIds;
use dbSNP::DBManager;
use Bio::EnsEMBL::Variation::Pipeline::VariantQC::RegisterDBSNPImport qw(register);

# If a config file was specified, parse it and override any specified options
my @opts;
if (scalar(@ARGV) == 1) {
  my $configfile = $ARGV[0];
  print STDOUT "Reading configuration file $configfile\n";
  open(CFG,'<',$configfile) or die("Could not open configuration file $configfile for reading");
  while (<CFG>) {
    chomp;
    next unless/\w+/;
    my ($name,$val) = (split/\s+/,$_,2);             #altered for cow which had space in primary assembly tag
    $val =~ s/\s+$//;
    push(@opts,('-' . $name,$val));
  }
  close(CFG);
  @ARGV = @opts;
}
  
my %options = ();
my @option_defs = (
  'species=s',
  'dbSNP_version=s',
  'shared_db=s',
  'tmpdir=s',
  'tmpfile=s',
  'limit=i',
  'mapping_file_dir=s',
  'schema_file=s',  ## ddl
  'schema_name:s',  ## postgreSQL schema name
  'dshost=s',
  'dsuser=s',
  'dspass=s',
  'dsport=i',
  'dsdbname=s',
  'registry_file=s',
  'mssql_driver=s',
  'skip_routine=s@',
  'scriptdir=s', 
  'logfile=s',
  'group_term=s',
  'group_label=s',
  'ensembl_version:s',
  'source_engine:s'
);

GetOptions(\%options,@option_defs);

debug("\n######### " . localtime() . " #########\n\tImport script launched\n");
print STDOUT "\n######### " . localtime() . " #########\n\tImport script launched\n";

my $LIMIT_SQL = $options{'limit'};
my $shared_db = $options{'shared_db'};
my $TMP_DIR = $options{'tmpdir'};
my $TMP_FILE = $options{'tmpfile'};
my $MAPPING_FILE_DIR = $options{'mapping_file_dir'};
my $SCHEMA_FILE = $options{'schema_file'};
my $schema_name = $options{'schema_name'};
my $GROUP_TERM  = $options{'group_term'};
my $GROUP_LABEL = $options{'group_label'};
my $species = $options{'species'};
my $dshost = $options{'dshost'};
my $dsuser = $options{'dsuser'};
my $dspass = $options{'dspass'};
my $dsport = $options{'dsport'};
my $dsdbname = $options{'dsdbname'};
my $registry_file = $options{'registry_file'};
my $mssql_driver = $options{'mssql_driver'};
my $scriptdir = $options{'scriptdir'};
my $logfile = $options{'logfile'};
my $ens_version = $options{'ensembl_version'};
my $dbSNP_version = $options{'dbSNP_version'};
my $source_engine = $options{'source_engine'} ||  'mysql';
my $ploidy = $options{'ploidy'} || 2;

my @skip_routines;
@skip_routines = @{$options{'skip_routine'}} if (defined($options{'skip_routine'}));

$ImportUtils::TMP_DIR = $TMP_DIR;
$ImportUtils::TMP_FILE = $TMP_FILE;

# Checking that some necessary arguments have been provided
die("Please provide a current schema definition file for the variation database, use -schema_file option!") unless (-e $SCHEMA_FILE);
# die("You must specify the dbSNP mirror host, port, user, pass, db and build version (-dshost, -dsport, -dsuser, -dspass, and -dbSNP_version options)") unless (defined($dshost) && defined($dsport) && defined($dsuser) && defined($dspass) && defined($dbSNP_BUILD_VERSION));
die("You must specify a temp dir and temp file (-tmpdir and -tmpfile options)") unless(defined($ImportUtils::TMP_DIR) && defined($ImportUtils::TMP_FILE));
die("You must specify the species. Use -species option") 
    unless (defined($species));
die("You must specify the dbSNP build. Use -dbSNP_version option") 
    unless (defined($dbSNP_version));
die("You must specify the dbSNP shared database. Use -shared_db option") 
    unless (defined($shared_db));
die("You must specify the sql driver, either through an environment variable (SYBASE) or the -sql_driver option") unless ((defined($mssql_driver) || defined($ENV{'SYBASE'})) ||  $source_engine =~/mysql|postgreSQL/i);

die("A schema name must be specified for postgreSQL imports")
    if  $source_engine =~/postgreSQL/i &&  !defined $schema_name;

warn("Note that the port for the dbSNP mirror is overridden by the freetds configuration file!\n") if (defined($dsport));


# Set the driver
$ENV{'SYBASE'} = $mssql_driver if (defined($mssql_driver));

# Set default option
$registry_file ||= $Bin . "/ensembl.registry";



# Open a file handle to the logfile. If it's not defined, use STDOUT
my $logh = *STDOUT;
if (defined($logfile)) {
  open(LOG,'>>',$logfile) or die ("Could not open logfile $logfile for writing");
  #Turn on autoflush for the logfile
  {
    my $ofh = select LOG;
    $| = 1;
    select $ofh;
  }
  $logh = *LOG;
  print $logh "\n######### " . localtime() . " #########\n\tImport script launched\n";
}

my $dbm = dbSNP::DBManager->new($registry_file,$species, $schema_name );
$dbm->dbSNP_shared($shared_db);
my $dbCore = $dbm->dbCore();

my ($cs) = @{$dbCore->get_CoordSystemAdaptor->fetch_all};
my $ASSEMBLY_VERSION = $cs->version();

#my $SPECIES_PREFIX = get_species_prefix($TAX_ID);

#create the dbSNP object for the specie we want to dump the data

my @parameters = (
  -DBManager => $dbm,
  -tmpdir => $TMP_DIR,
  -tmpfile => $TMP_FILE,
  -limit => $LIMIT_SQL,
  -mapping_file_dir => $MAPPING_FILE_DIR,
  -dbSNP_version => $dbSNP_version,
  -assembly_version => $ASSEMBLY_VERSION,
  -group_term  => $GROUP_TERM,
  -group_label => $GROUP_LABEL,
  -skip_routines => \@skip_routines,
  -scriptdir => $scriptdir,
  -log => $logh,
  -source_engine => $source_engine,
  -schema_name   => $schema_name
);

my $import_object;

if ($species =~ m/felix_cattus/i ||  $species =~ m/tetraodon/i) {
  $import_object = dbSNP::MappingChromosome->new(@parameters);
}
elsif ($species =~ m/zebrafish|danio/i || 
       $species =~ m/chimp|troglodytes/i || 
       $species =~ m/dog|canis/i     ||
       $species =~ m/gallus_gallus/i || 
       $species =~ m/rat/i ||
       $species =~ m/horse|equus/i ||  
       $species =~ m/opossum/i ||
       $species =~ m/mus_musculus/i || 
       $species =~ m/bos_taurus/i   || 
       $species =~ m/sus_scrofa/i   ||  
       $species =~ m/felix_cattus/i || 
       $species =~ m/zebrafinch|taeniopygia_guttata/i || 
       $species =~ m/tetraodon/i || 
       $species =~ m/orangutan|Pongo_abelii/i || 
       $species =~ m/monodelphis_domestica/i  || 
       $species =~ m/macaca_mulatta/i   ||
       $species =~ m/ovis_aries/i       ||
       $species =~ m/zebrafinch|taeniopygia_guttata/i ||
       $species =~ m/cat|felis_catus/
    ) {
    $import_object = dbSNP::GenericChromosome->new(@parameters);
}
elsif ($species =~ m/platypus|anatinus/i) {
  $import_object = dbSNP::GenericContig->new(@parameters);
}
elsif ($species =~ m/mosquitos/i) {
  $import_object = dbSNP::Mosquito->new(@parameters);
}
elsif ($species =~ m/human|homo/i) {
  $import_object = dbSNP::EnsemblIds->new(@parameters);
}
else {
  die("The species needs to have a module to use for import hardcoded into this script");
}

$import_object->{'schema_file'} = $SCHEMA_FILE;

my $clock = Progress->new();
$clock->checkpoint('start_dump');

$import_object->dump_dbSNP();

$clock->checkpoint('end_dump');
print $logh $clock->duration('start_dump','end_dump');

### Previously this script copied tmp_individual_genotypes to alleles for mouse, rat, chicken and dog.
### This behaviour ceased as of 30/1/2013



##d update meta 
my $meta_ins_sth = $dbm->dbVar()->dbc->db_handle->prepare(qq[ INSERT ignore INTO meta (species_id, meta_key, meta_value) values (?,?,?)]);
my $meta_upd_sth = $dbm->dbVar()->dbc->db_handle->prepare(qq[ UPDATE meta set meta_value = ? where  meta_key = ?]);

$meta_ins_sth->execute('1', 'species.production_name', $dbm->dbVar()->species() ) ||die;
$meta_ins_sth->execute('1', 'ploidy', $ploidy ) ||die;


if (defined $ens_version){
    $meta_upd_sth->execute($ens_version, 'schema_version' ) ||die;
}


### update production db as final step

my $dbSNP_name = $dbm->dbSNP()->dbc->dbname();
$dbSNP_name =~ s/\_\d+$//;

my %data;

$data{dbSNP_name}       = $dbSNP_name;
$data{dbSNP_name} =~ s/^dbSNP\_//; # prefix used for mysql mirror
$data{species}          = $dbm->dbVar()->species();
$data{ensdb_name}       = $dbm->dbVar()->dbc->dbname();
$data{registry}         = $dbm->registry();
$data{ensdb_version}    = $ens_version;
$data{assembly_version} = $ASSEMBLY_VERSION;
$data{pwd}              = $TMP_DIR;
register(\%data);






debug(localtime() . "\tAll done!");

#Close the filehandle to the logfile if one was specified
close($logh) if (defined($logfile));
  

sub usage {
    my $msg = shift;
    
    print STDERR <<EOF;
    
  usage: perl dbSNP.pl <options>
      
    options:
      -dshost <hostname>          hostname of dbSNP MSSQL database (default = dbsnp)
      -dsuser <user>              username of dbSNP MSSQL database (default = dbsnpro)
      -dspass <pass>              password of dbSNP MSSQL database
      -dsport <port>              TCP port of dbSNP MSSQL database (default = 1026)
      -dsdbname <dbname>          dbname of dbSNP MySQL database   (default = dbSNP_121)
      -chost <hostname>           hostname of core Ensembl MySQL database (default = ecs2)
      -cuser <user>               username of core Ensembl MySQL database (default = ensro)
      -cpass <pass>               password of core Ensembl MySQL database
      -cport <port>               TCP port of core Ensembl MySQL database (default = 3364)
      -cdbname <dbname>           dbname of core Ensembl MySQL database
      -vhost <hostname>           hostname of variation MySQL database to write to
      -vuser <user>               username of variation MySQL database to write to (default = ensadmin)
      -vpass <pass>               password of variation MySQL database to write to
      -vport <port>               TCP port of variation MySQL database to write to (default = 3306)
      -vdbname <dbname>           dbname of variation MySQL database to write to
      -limit <number>             limit the number of rows transfered for testing
      -tmpdir <dir>               temporary directory to use (with lots of space!)
      -tmpfile <filename>         temporary filename to use
      -mapping_file <filename>    file containing the mapping data
      -group_term <group_term>    select the group_term to import
      -group_label <group_label>  select the group_label to import
                        
EOF

      die("\n$msg\n\n");
}

