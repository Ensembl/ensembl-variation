#!/usr/local/ensembl/bin/perl -w
# this is the script to fill the new variation 
# schema with data from dbSNP
# we use the local mysql copy of dbSNP at the sanger center
# the script will call the dbSNP factory to create the object that will deal with
# the creation of the data according to the specie

use strict;
use warnings;
use DBI;
use DBH;
use Getopt::Long;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use ImportUtils qw(dumpSQL debug create_and_load load);
use Bio::EnsEMBL::Utils::Exception qw(throw);


use dbSNP::GenericContig;
use dbSNP::GenericChromosome;
use dbSNP::Zebrafish;
use dbSNP::Mosquito;
use dbSNP::Human;

my ($TAX_ID, $LIMIT_SQL, $CONTIG_SQL, $TMP_DIR, $TMP_FILE, $ALLDIFF_FILE);

my $dbSNP;
my $dbVar;
my $dbCore;

{
  my($dshost, $dsuser, $dspass, $dsport, $dsdbname, # dbSNP db
     $chost, $cuser, $cpass, $cport, $cdbname,      # ensembl core db
     $vhost, $vuser, $vpass, $vport, $vdbname,      # ensembl variation db
     $limit);

  GetOptions('dshost=s'   => \$dshost,
             'dsuser=s'   => \$dsuser,
             'dspass=s'   => \$dspass,
             'dsport=i'   => \$dsport,
             'dsdbname=s' => \$dsdbname,
             'chost=s'   => \$chost,
             'cuser=s'   => \$cuser,
             'cpass=s'   => \$cpass,
             'cport=i'   => \$cport,
             'cdbname=s' => \$cdbname,
             'vhost=s'   => \$vhost,
             'vuser=s'   => \$vuser,
             'vpass=s'   => \$vpass,
             'vport=i'   => \$vport,
             'vdbname=s' => \$vdbname,
             'tmpdir=s'  => \$ImportUtils::TMP_DIR,
             'tmpfile=s' => \$ImportUtils::TMP_FILE,
             'limit=i'   => \$limit,
	     'alldiff=s' => \$ALLDIFF_FILE);

  $dshost   ||= 'cbi2.internal.sanger.ac.uk';
  $dsdbname ||= 'dbSNP_123';
  $dsuser   ||= 'dbsnpro';
  $dsport   ||= 3306;

  $chost    ||= 'ecs2';
  $cuser    ||= 'ensro';
  $cport    ||= 3365;

  $vhost    ||= 'ecs2';
  $vport    ||= 3362;
  $vuser    ||= 'ensadmin';

  usage('-cdbname argument is required.') if(!$cdbname);
  usage('-vdbname argument is required.') if(!$vdbname);

  $TMP_DIR  = $ImportUtils::TMP_DIR;
  $TMP_FILE = $ImportUtils::TMP_FILE;


  $LIMIT_SQL = ($limit) ? " LIMIT $limit " : '';

  $dbSNP = DBH->connect
    ("DBI:mysql:host=$dshost;dbname=$dsdbname;port=$dsport",$dsuser, $dspass,
    {'RaiseError' => 1});
  die("Could not connect to dbSNP db: $!") if(!$dbSNP);

  $dbSNP->{mysql_auto_reconnect} = 1;

  $dbSNP->do("SET SESSION wait_timeout = 2678200");

  $dbVar = DBH->connect
      ("DBI:mysql:host=$vhost;dbname=$vdbname;port=$vport",$vuser, $vpass,
       {'RaiseError' => 1});
  die("Could not connect to variation database: $!") if(!$dbVar);

  $dbCore = Bio::EnsEMBL::DBSQL::DBAdaptor->new
    (-host   => $chost,
     -user   => $cuser,
     -pass   => $cpass,
     -port   => $cport,
     -dbname => $cdbname);

  my $mc = $dbCore->get_MetaContainer();
  my $species = $mc->get_Species();

  throw("Unable to determine species from core database.") if(!$species);

  $TAX_ID = $mc->get_taxonomy_id();


  if($species->binomial() eq 'Homo sapiens') {
#    $CONTIG_SQL = ' CONCAT(contig_acc, ".", contig_ver) ';
    $CONTIG_SQL = ' contig_acc ';
  } else {
    $CONTIG_SQL = ' contig_acc ';
  }
}

my $SPECIES_PREFIX = get_species_prefix($TAX_ID);

#create the dbSNP object for the specie we want to dump the data

if ($SPECIES_PREFIX eq 'dr'){
    if (!$ALLDIFF_FILE) {throw ("file with mappings not provided")}
    #danio rerio (zebra-fish)
    my $zebrafish = dbSNP::Zebrafish->new(-dbSNP => $dbSNP,
					  -dbCore => $dbCore,
					  -dbVariation => $dbVar,					   
					  -tmpdir => $TMP_DIR,
					  -tmpfile => $TMP_FILE,
					  -limit => $LIMIT_SQL,
					  -alldiff => $ALLDIFF_FILE,
					  -taxID => $TAX_ID
					  );
    $zebrafish->dump_dbSNP();
}
elsif ($SPECIES_PREFIX eq 'mm'){
    #mus-musculus (mouse)
    my $mouse = dbSNP::GenericContig->new(-dbSNP => $dbSNP,
					  -dbCore => $dbCore,
					  -dbVariation => $dbVar,
					  -tmpdir => $TMP_DIR,
					  -tmpfile => $TMP_FILE,
					  -limit => $LIMIT_SQL,
					  -taxID => $TAX_ID,
					  -contigSQL => $CONTIG_SQL,
					  -species_prefix => $SPECIES_PREFIX
					  );
    $mouse->dump_dbSNP();
}
elsif ($SPECIES_PREFIX eq 'gga'){
    #gallus gallus (chicken)
    my $chicken = dbSNP::GenericChromosome->new(-dbSNP => $dbSNP,
						-dbCore => $dbCore,
						-dbVariation => $dbVar,
						-tmpdir => $TMP_DIR,
						-tmpfile => $TMP_FILE,
						-limit => $LIMIT_SQL,
						-taxID => $TAX_ID,
						-species_prefix => $SPECIES_PREFIX
					      );
    $chicken->dump_dbSNP();
}
elsif ($SPECIES_PREFIX eq 'rn'){
   #rattus norvegiccus (rat)
    my $rat = dbSNP::GenericChromosome->new(-dbSNP => $dbSNP,
					    -dbCore => $dbCore,
					    -dbVariation => $dbVar,
					    -tmpdir => $TMP_DIR,
					    -tmpfile => $TMP_FILE,
					    -limit => $LIMIT_SQL,
					    -taxID => $TAX_ID,
					    -species_prefix => $SPECIES_PREFIX
					    );
    $rat->dump_dbSNP();
}
elsif ($SPECIES_PREFIX eq 'ag'){
    #(mosquito)
    my $mosquito = dbSNP::Mosquito->new(-dbSNP => $dbSNP,
					-dbCore => $dbCore,
					-dbVariation => $dbVar,
					-tmpdir => $TMP_DIR,
					-tmpfile => $TMP_FILE,
					-limit => $LIMIT_SQL,
					-taxID => $TAX_ID,
					-species_prefix => $SPECIES_PREFIX);
    $mosquito->dump_dbSNP();
}
else{
    #homo sa/piens (human)
    my $human = dbSNP::Human->new(-dbSNP => $dbSNP,
					      -dbCore => $dbCore,
					      -dbVariation => $dbVar,
					      -tmpdir => $TMP_DIR,
					      -tmpfile => $TMP_FILE,
					      -limit => $LIMIT_SQL,
					      -taxID => $TAX_ID,
					      -contigSQL => $CONTIG_SQL,
					      -species_prefix => $SPECIES_PREFIX
					      );
    $human->dump_dbSNP();
}


sub get_species_prefix {
  my $tax_id = shift;

  my $arr = $dbSNP->selectall_arrayref
    (qq{SELECT ot.prefix
        FROM   OrganismTax ot
        WHERE  ot.tax_id = $tax_id});

  if(@$arr) {
    return $arr->[0]->[0];
  }

  warn("tax_id=$tax_id not found in OrganismTax table." .
       "Assuming no species prefix");
  return '';
}


sub usage {
    my $msg = shift;
    
    print STDERR <<EOF;
    
  usage: perl dbSNP.pl <options>
      
    options:
      -dshost <hostname>   hostname of dbSNP MySQL database (default = cbi2.internal.sanger.ac.uk)
      -dsuser <user>       username of dbSNP MySQL database (default = dbsnpro)
      -dspass <pass>       password of dbSNP MySQL database
      -dsport <port>       TCP port of dbSNP MySQL database (default = 3306)
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
      -alldiff <filename>  file containing the mapping data for the zebrafish specie
EOF

      die("\n$msg\n\n");
}

