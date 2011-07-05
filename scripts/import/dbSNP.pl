#!/usr/local/bin/perl -w
# this is the script to fill the new variation 
# schema with data from dbSNP
# we use the local mssql copy of dbSNP at the sanger center
# the script will call the dbSNP factory to create the object that will deal with
# the creation of the data according to the specie

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

# If a config file was specified, parse it and override any specified options
my @opts;
if (scalar(@ARGV) == 1) {
  my $configfile = $ARGV[0];
  print STDOUT "Reading configuration file $configfile\n";
  open(CFG,'<',$configfile) or die("Could not open configuration file $configfile for reading");
  while (<CFG>) {
    chomp;
    my ($name,$val) = $_ =~ m/^\s*(\S+)\s+([^\s]+)/;
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
  'schema_file=s',
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
	'group_label=s'
);

GetOptions(\%options,@option_defs);

debug("\n######### " . localtime() . " #########\n\tImport script launched\n");
print STDOUT "\n######### " . localtime() . " #########\n\tImport script launched\n";

my $LIMIT_SQL = $options{'limit'};
my $dbSNP_BUILD_VERSION = $options{'dbSNP_version'};
my $shared_db = $options{'shared_db'};
my $TMP_DIR = $options{'tmpdir'};
my $TMP_FILE = $options{'tmpfile'};
my $MAPPING_FILE_DIR = $options{'mapping_file_dir'};
my $SCHEMA_FILE = $options{'schema_file'};
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
my @skip_routines;
@skip_routines = @{$options{'skip_routine'}} if (defined($options{'skip_routine'}));

$ImportUtils::TMP_DIR = $TMP_DIR;
$ImportUtils::TMP_FILE = $TMP_FILE;

# Checking that some necessary arguments have been provided
die("Please provide a current schema definition file for the variation database, use -schema_file option!") unless (-e $SCHEMA_FILE);
# die("You must specify the dbSNP mirror host, port, user, pass, db and build version (-dshost, -dsport, -dsuser, -dspass, and -dbSNP_version options)") unless (defined($dshost) && defined($dsport) && defined($dsuser) && defined($dspass) && defined($dbSNP_BUILD_VERSION));
die("You must specify a temp dir and temp file (-tmpdir and -tmpfile options)") unless(defined($ImportUtils::TMP_DIR) && defined($ImportUtils::TMP_FILE));
die("You must specify the species. Use -species option") unless (defined($species));
die("You must specify the dbSNP build. Use -dbSNP_version option") unless (defined($dbSNP_BUILD_VERSION));
die("You must specify the dbSNP shared database. Use -shared_db option") unless (defined($shared_db));
die("You must specify the sql driver, either through an environment variable (SYBASE) or the -sql_driver option") unless (defined($mssql_driver) || defined($ENV{'SYBASE'}));

warn("Note that the port for the dbSNP mirror is overridden by the freetds configuration file!\n") if (defined($dsport));
warn("Make sure you have a updated ensembl.registry file!\n");

# Set the driver
$ENV{'SYBASE'} = $mssql_driver if (defined($mssql_driver));

# Set default option
$registry_file ||= $Bin . "/ensembl.registry";

# Open a file handle to the logfile. If it's not defined, use STDOUT
my $logh = *STDOUT;
if (defined($logfile)) {
  open(LOG,'>>',$logfile) or die ("Could not open logfile $logfile for writing");
  #�Turn on autoflush for the logfile
  {
    my $ofh = select LOG;
    $| = 1;
    select $ofh;
  }
  $logh = *LOG;
  print $logh "\n######### " . localtime() . " #########\n\tImport script launched\n";
}

=head
Bio::EnsEMBL::Registry->load_all( $registry_file );

my $cdba = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'core') or die ("Could not get core DBadaptor");
my $vdba = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'variation') or die ("Could not get DBadaptor to destination variation database");
my $snpdba = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'dbsnp') or die ("Could not get DBadaptor to dbSNP source database");
# Set the disconnect_when_inactive property
#$cdba->dbc->disconnect_when_inactive(1);
#$vdba->dbc->disconnect_when_inactive(1);
#$snpdba->dbc->disconnect_when_inactive(1);

$vdba->dbc->{mysql_auto_reconnect} = 1;
$cdba->dbc->{mysql_auto_reconnect} = 1;

$vdba->dbc->do("SET SESSION wait_timeout = 2678200");
$cdba->dbc->do("SET SESSION wait_timeout = 2678200");

 Set some variables on the MySQL server that can speed up table read/write/loads
my $stmt = qq{
  SET SESSION
    bulk_insert_buffer_size=512*1024*1024
};
$vdba->dbc->do($stmt);

if (!$dsdbname) {
  my $TAX_ID = $cdba->get_MetaContainer()->get_taxonomy_id() or throw("Unable to determine taxonomy id from core database for species $species.");
  my $version = $dbSNP_BUILD_VERSION;
  $version =~ s/^b//;
 $dsdbname = "$species\_$TAX_ID\_$version";
}

my $dbSNP = $snpdba->dbc;
my $dbVar = $vdba->dbc;
my $dbCore = $cdba;

#my $my_species = $mc->get_Species();
=cut

my $dbm = dbSNP::DBManager->new($registry_file,$species);
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
  -dbSNP_version => $dbSNP_BUILD_VERSION,
  -assembly_version => $ASSEMBLY_VERSION,
  -group_term  => $GROUP_TERM,
	-group_label => $GROUP_LABEL,
  -skip_routines => \@skip_routines,
  -scriptdir => $scriptdir,
  -log => $logh
);
my $import_object;
if ($species =~ m/zebrafish/i || $species =~ m/^pig$/i || $species eq 'cat' || $species =~ m/zebrafinch/i || $species =~ m/tetraodon/i) {
  $import_object = dbSNP::MappingChromosome->new(@parameters);
}
elsif ($species =~ m/chimp/i || $species =~ m/chicken/i || $species =~ m/rat/i || $species =~ m/horse/i || $species =~ m/platypus/i || $species =~ m/opossum/i || $species =~ m/mouse/i || $species =~ m/cow/i) {
  $import_object = dbSNP::GenericChromosome->new(@parameters);
}
elsif ($species =~ m/dog/i) {
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

if ($species =~ m/mouse/i || $species =~ m/chicken/i || $species =~ m/rat/i || $species =~ m/dog/i) {
  # Add strain information
  $clock->checkpoint('add_strains');
  
  add_strains($dbm->dbVar());
  
  print $logh $clock->duration();
}

debug(localtime() . "\tAll done!");

#�Close the filehandle to the logfile if one was specified
close($logh) if (defined($logfile));
  

#gets all the individuals and copies the data in the Population table as strain, and the individual_genotype as allele
sub add_strains{
  my $dbVariation = shift;

  debug(localtime() . "\tIn add_strains...");

  #first, copy the data from Individual to Population, necessary to group by name to remove duplicates
  # (some individuals have the same name, but different individual_id)
  $dbVariation->do(qq{ALTER TABLE sample ADD COLUMN is_strain int});

  $dbVariation->do(qq{INSERT INTO sample (name, description,is_strain)
		      SELECT s.name, s.description,1
		      FROM individual i, sample s
		      WHERE i.sample_id = s.sample_id
		      GROUP BY s.name
		     });
  #and insert the new strain in the population table
  $dbVariation->do(qq{INSERT INTO population (sample_id, is_strain)
		      SELECT sample_id, is_strain
		      FROM sample
		      WHERE is_strain = 1
		     });
  #then, populate the individual_population table with the relation between individual and strains
  $dbVariation->do(qq{INSERT INTO individual_population (individual_sample_id, population_sample_id)
		      SELECT s1.sample_id, s2.sample_id
		      FROM sample s1, sample s2, individual i
		      WHERE s1.name = s2.name
		      AND s2.is_strain = 1
		      AND s1.sample_id = i.sample_id
		     });

  #check to see allele_1 and allele_2 is same or not to determine frequency

  for my $gtype_table ('tmp_individual_genotype_single_bp','individual_genotype_multiple_bp') {
    debug(localtime() . "\tcreate tmp_freq table...from $gtype_table");

    $dbVariation->do(qq{CREATE TABLE tmp_freq (variation_id int, allele_1 text, sample_id int, gtype char (4), key (sample_id))});
    #add allele_1 into tmp_freq
    $dbVariation->do(qq{INSERT INTO tmp_freq
			SELECT distinct variation_id, allele_1, sample_id,
			IF (allele_1 != allele_2, 'DIFF', 'SAME')
			FROM $gtype_table
		       });
    #add allele_2 into tmp_freq only if it's different from allele_1
    $dbVariation->do(qq{INSERT INTO tmp_freq
			SELECT distinct variation_id, allele_2, sample_id, 'DIFF'
			FROM $gtype_table
			WHERE allele_1 != allele_2
		       });

    #and finally, add the alleles to the Allele table for the strains
    debug(localtime() . "\tInsert into allele in add_strain...");
    $dbVariation->do(qq{INSERT INTO allele (variation_id, sample_id, allele, frequency)
			SELECT tf.variation_id, ip.population_sample_id, tf.allele_1,
			IF (tf.gtype = 'SAME',1,0.5)
			FROM tmp_freq tf, individual_population ip, sample s
			WHERE tf.sample_id = ip.individual_sample_id
			AND   ip.population_sample_id = s.sample_id
			AND   s.is_strain = 1
		       });

    $dbVariation->do(qq{DROP TABLE tmp_freq});
  }

  $dbVariation->do(qq{ALTER TABLE sample DROP COLUMN is_strain});
}


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

