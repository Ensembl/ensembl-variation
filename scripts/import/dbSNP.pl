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

my ($TAX_ID, $LIMIT_SQL, $dbSNP_BUILD_VERSION, $TMP_DIR, $TMP_FILE, $MAPPING_FILE_DIR, $MASTER_SCHEMA_DB);

my ($species,$limit);
my($dshost, $dsuser, $dspass, $dsport, $dsdbname);
my $registry_file;
my $sql_driver;
my @skip_routines;
my $logfile;

GetOptions('species=s'      => \$species,
	   'dbSNP_version=s'=> \$dbSNP_BUILD_VERSION, ##such as b125
	   'tmpdir=s'       => \$ImportUtils::TMP_DIR,
	   'tmpfile=s'      => \$ImportUtils::TMP_FILE,
	   'limit=i'        => \$limit,
	   'mapping_file_dir=s' => \$MAPPING_FILE_DIR,
	   'master_schema_db=s' => \$MASTER_SCHEMA_DB,
	   'dshost=s' => \$dshost,
	   'dsuser=s' => \$dsuser,
	   'dspass=s' => \$dspass,
	   'dsport=i' => \$dsport,
	   'dsdbname=s' => \$dsdbname,
	   'registry_file=s' => \$registry_file,
	   'sql_driver=s' => \$sql_driver,
	   'skip_routine=s' => \@skip_routines,
	   'logfile=s' => \$logfile
	  );

debug("\n######### " . localtime() . " #########\n\tImport script launched\n");
print STDOUT "\n######### " . localtime() . " #########\n\tImport script launched\n";

# Checking that some necessary arguments have been provided
die("Must know the master schema database, use -master_schema_db option!") unless (defined($MASTER_SCHEMA_DB));
die("You must specify the dbSNP mirror host, user, pass, db and build version (-dshost, -dsuser, -dspass, and -dbSNP_version options)") unless (defined($dshost) && defined($dsuser) && defined($dspass) && defined($dbSNP_BUILD_VERSION));
die("You must specify a temp dir and temp file (-tmpdir and -tmpfile options)") unless(defined($ImportUtils::TMP_DIR) && defined($ImportUtils::TMP_FILE));
die("You must specify the species. Use -species option") unless (defined($species));
die("You must specify the sql driver, either through an environment variable (SYBASE) or the -sql_driver option") unless (defined($sql_driver) || defined($ENV{'SYBASE'}));

warn("Note that the port for the dbSNP mirror is overridden by the freetds configuration file!\n") if (defined($dsport));
warn("Make sure you have a updated ensembl.registry file!\n");

# Set the driver
$ENV{'SYBASE'} = $sql_driver if (defined($sql_driver));

# Set default option
$registry_file ||= $Bin . "/ensembl.registry";

# Open a file handle to the logfile. If it's not defined, use STDOUT
my $logh = *STDOUT;
if (defined($logfile)) {
  open(LOG,'>>',$logfile) or die ("Could not open logfile $logfile for writing");
  #ÊTurn on autoflush for the logfile
  {
    my $ofh = select LOG;
    $| = 1;
    select $ofh;
  }
  $logh = *LOG;
  print $logh "\n######### " . localtime() . " #########\n\tImport script launched\n";
}

Bio::EnsEMBL::Registry->load_all( $registry_file );

my $cdba = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'core') or die ("Could not get core DBadaptor");
my $vdba = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'variation') or die ("Could not get DBadaptor to destination variation database");
$vdba->dbc->{mysql_auto_reconnect} = 1;
$vdba->dbc->do("SET SESSION wait_timeout = 2678200");
#my $sdba = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'dbsnp');

$TAX_ID = $cdba->get_MetaContainer()->get_taxonomy_id();
throw("Unable to determine taxonomy id from core database for species $species.") unless $TAX_ID;

if (!$dsdbname) {
  my $version = $dbSNP_BUILD_VERSION;
  $version =~ s/^b//;
  $dsdbname = "$species\_$TAX_ID\_$version";
}

my $source = 'DBI:Sybase:';

# Note that the port is specified under [dbsnp] in the freetds.conf file!
# Server should be 'dbsnp' which then points the database driver to the corresponding entry in the freetds.conf file.
# It can be overridden with a local .freetds.conf file in the user's home directory
my $dbSNP = DBI->connect($source . 
                         "timeout=2678200;".
                         "database=$dsdbname;".
                         "server=$dshost",
                          $dsuser,
                          $dspass
                        ) or die("Could not connect to database $DBI::errstr");

my $dbVar = $vdba->dbc;
my $dbCore = $cdba;

###find a large tmp directory to contain the tmp file
$TMP_DIR  = $ImportUtils::TMP_DIR;
$TMP_FILE = $ImportUtils::TMP_FILE;

$LIMIT_SQL = $limit;


#my $my_species = $mc->get_Species();

my ($cs) = @{$dbCore->get_CoordSystemAdaptor->fetch_all};
my $ASSEMBLY_VERSION = $cs->version();

#my $SPECIES_PREFIX = get_species_prefix($TAX_ID);

#create the dbSNP object for the specie we want to dump the data

my @parameters = (
  -dbSNP => $dbSNP,
  -dbCore => $dbCore,
  -dbVar => $dbVar,
  -snp_dbname => $dsdbname,
  -tmpdir => $TMP_DIR,
  -tmpfile => $TMP_FILE,
  -limit => $LIMIT_SQL,
  -mapping_file_dir => $MAPPING_FILE_DIR,
  -dbSNP_version => $dbSNP_BUILD_VERSION,
  -assembly_version => $ASSEMBLY_VERSION,
  -skip_routines => \@skip_routines,
  -log => $logh
);
my $import_object;
if ($species =~ m/fish/i || $species =~ m/^pig$/i) {
  $import_object = dbSNP::MappingChromosome->new(@parameters);
}
elsif ($species =~ m/chimp/i || $species =~ m/chicken/i || $species =~ m/rat/i || $species =~ m/horse/i || $species =~ m/platypus/i) {
  $import_object = dbSNP::GenericChromosome->new(@parameters);
}
elsif ($species =~ m/mouse/i || $species =~ m/dog/i) {
  $import_object = dbSNP::GenericContig->new(@parameters);
}
elsif ($species =~ m/mosquitos/i) {
  $import_object = dbSNP::Mosquito->new(@parameters);
}
elsif ($species =~ m/human|homo/i) {
  $import_object = dbSNP::EnsemblIds->new(@parameters);
}
$import_object->{'master_schema_db'} = $MASTER_SCHEMA_DB;

my $clock = Progress->new();
$clock->checkpoint('start_dump');

$import_object->dump_dbSNP();

$clock->checkpoint('end_dump');
print $logh $clock->duration('start_dump','end_dump');

if ($species =~ m/mouse/i || $species =~ m/chicken/i || $species =~ m/rat/i || $species =~ m/dog/i) {
  # Add strain information
  $clock->checkpoint('add_strains');
  
  add_strains($dbVar);
  
  print $logh $clock->duration();
}

debug(localtime() . "\tAll done!");

#ÊClose the filehandle to the logfile if one was specified
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
      -dshost <hostname>   hostname of dbSNP MSSQL database (default = dbsnp)
      -dsuser <user>       username of dbSNP MSSQL database (default = dbsnpro)
      -dspass <pass>       password of dbSNP MSSQL database
      -dsport <port>       TCP port of dbSNP MSSQL database (default = 1026)
      -dsdbname <dbname>   dbname of dbSNP MySQL database   (default = dbSNP_121)
      -chost <hostname>    hostname of core Ensembl MySQL database (default = ecs2)
      -cuser <user>        username of core Ensembl MySQL database (default = ensro)
      -cpass <pass>        password of core Ensembl MySQL database
      -cport <port>        TCP port of core Ensembl MySQL database (default = 3364)
      -cdbname <dbname>    dbname of core Ensembl MySQL database
      -vhost <hostname>    hostname of variation MySQL database to write to
      -vuser <user>        username of variation MySQL database to write to (default = ensadmin)
      -vpass <pass>        password of variation MySQL database to write to
      -vport <port>        TCP port of variation MySQL database to write to (default = 3306)
      -vdbname <dbname>    dbname of variation MySQL database to write to
      -limit <number>      limit the number of rows transfered for testing
      -tmpdir <dir>        temporary directory to use (with lots of space!)
      -tmpfile <filename>  temporary filename to use
      -mapping_file <filename>  file containing the mapping data
EOF

      die("\n$msg\n\n");
}

