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

my ($TAX_ID, $LIMIT_SQL, $TMP_DIR, $TMP_FILE, $ALLDIFF_FILE);

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
  $dsdbname ||= 'dbSNP_124';
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
					  -species_prefix => $SPECIES_PREFIX
					  );
    $mouse->dump_dbSNP();
    add_strains($dbVar); #add strain information
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
    add_strains($dbVar); #add strain information
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
    add_strains($dbVar); #add strain information
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
elsif ($SPECIES_PREFIX eq 'cfa'){
    #cannis familiaris (dog)
    my $dog = dbSNP::GenericContig->new(-dbSNP => $dbSNP,
					  -dbCore => $dbCore,
					  -dbVariation => $dbVar,
					  -tmpdir => $TMP_DIR,
					  -tmpfile => $TMP_FILE,
					  -limit => $LIMIT_SQL,
					  -taxID => $TAX_ID,
					  -species_prefix => $SPECIES_PREFIX
					  );
    $dog->dump_dbSNP();
    add_strains($dbVar); #add strain information

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

#gets all the individuals and copies the data in the Population table as strain, and the individual_genotype as allele
sub add_strains{
    my $dbVariation = shift;

    #first, copy the data from Individual to Population, necessary to group by name to remove duplicates
    # (some individuals have the same name, but different individual_id)
    $dbVariation->do(qq{INSERT INTO population (name, description, is_strain)
				      SELECT name, description, 1
				      FROM individual
				      GROUP BY name
				  });

    #then, populate the individual_population table with the relation between individual and strains
    $dbVariation->do(qq{INSERT INTO individual_population (individual_id, population_id)
				      SELECT i.individual_id, p.population_id
				      FROM individual i, population p
				      WHERE i.name = p.name
				      AND p.is_strain = 1
				  });

    #and finally, add the alleles to the Allele table for the strains
    #first, insert the single_bp alleles
    $dbVariation->do(qq{INSERT INTO allele (variation_id, population_id, allele, frequency)
				      SELECT ig.variation_id, ip.population_id, ig.allele_1, 1
				      FROM individual_genotype_single_bp ig, individual_population ip, population p
				      WHERE ig.individual_id = ip.individual_id
				      AND   ip.population_id = p.population_id
				      AND   p.is_strain = 1
				  });
    #do the same for the multiple_bp alleles
    $dbVariation->do(qq{INSERT INTO allele (variation_id, population_id, allele, frequency)
				      SELECT ig.variation_id, ip.population_id, ig.allele_1, 1
				      FROM individual_genotype_multiple_bp ig, individual_population ip, population p
				      WHERE ig.individual_id = ip.individual_id
				      AND   ip.population_id = p.population_id
				      AND   p.is_strain = 1
				  });
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

