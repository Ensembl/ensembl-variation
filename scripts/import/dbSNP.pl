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
use Bio::EnsEMBL::Registry;
use FindBin qw( $Bin );

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use ImportUtils qw(dumpSQL debug create_and_load load);
use Bio::EnsEMBL::Utils::Exception qw(throw);


use dbSNP::GenericContig;
use dbSNP::GenericChromosome;
use dbSNP::MappingChromosome;
use dbSNP::Mosquito;
use dbSNP::Human;

my ($TAX_ID, $LIMIT_SQL, $TMP_DIR, $TMP_FILE, $MAPPING_FILE);

my ($species,$limit);
my($dshost, $dsuser, $dspass, $dsport, $dsdbname);

GetOptions('species=s'      => \$species,
	   'tmpdir=s'       => \$ImportUtils::TMP_DIR,
	   'tmpfile=s'      => \$ImportUtils::TMP_FILE,
	   'limit=i'        => \$limit,
	   'mapping_file=s' => \$MAPPING_FILE);

warn("Make sure you have a updated ensembl.registry file!\n");

my $registry_file ||= $Bin . "/ensembl.registry";

Bio::EnsEMBL::Registry->load_all( $registry_file );

my $cdba = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'core');
my $vdba = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'variation');
my $sdba = Bio::EnsEMBL::Registry->get_DBAdaptor('dbsnp','dbsnp');


my $dbSNP = $sdba->dbc->db_handle; # Equivalent to direct connection
my $dbVar = $vdba->dbc->db_handle;
my $dbCore = $cdba;

$TMP_DIR  = $ImportUtils::TMP_DIR;
$TMP_FILE = $ImportUtils::TMP_FILE;


$LIMIT_SQL = ($limit) ? " LIMIT $limit " : '';

$dbSNP->{mysql_auto_reconnect} = 1;

$dbSNP->do("SET SESSION wait_timeout = 2678200");

my $mc = $cdba->get_MetaContainer();
my $my_species = $mc->get_Species();

throw("Unable to determine species from core database.") if(!$my_species);

$TAX_ID = $mc->get_taxonomy_id();




my $SPECIES_PREFIX = get_species_prefix($TAX_ID);

#create the dbSNP object for the specie we want to dump the data
if ($SPECIES_PREFIX eq 'dr'){
  if (!$MAPPING_FILE) {throw ("file with mappings not provided")}
  #danio rerio (zebra-fish)
  my $zebrafish = dbSNP::MappingChromosome->new(-dbSNP => $dbSNP,
						-dbCore => $dbCore,
						-dbVariation => $dbVar,					                                          -snp_dbname => $sdba->dbc->dbname,
						-tmpdir => $TMP_DIR,
						-tmpfile => $TMP_FILE,
						-limit => $LIMIT_SQL,
						-mapping_file => $MAPPING_FILE,
						-taxID => $TAX_ID
					       );
  $zebrafish->dump_dbSNP();
}
elsif ($SPECIES_PREFIX eq 'mm'){
  #mus-musculus (mouse)
  #my $mouse = dbSNP::GenericContig->new(-dbSNP => $dbSNP,
  my $mouse = dbSNP::MappingChromosome->new(-dbSNP => $dbSNP,
					    -dbCore => $dbCore,
					    -dbVariation => $dbVar,
					    -snp_dbname => $sdba->dbc->dbname,
					    -tmpdir => $TMP_DIR,
					    -tmpfile => $TMP_FILE,
					    -limit => $LIMIT_SQL,
					    -taxID => $TAX_ID,
					    -species_prefix => $SPECIES_PREFIX,
					    -mapping_file => $MAPPING_FILE,
					   );
  $mouse->dump_dbSNP();
  add_strains($dbVar); #add strain information
}
elsif ($SPECIES_PREFIX eq 'gga'){
  #gallus gallus (chicken)
  my $chicken = dbSNP::GenericChromosome->new(-dbSNP => $dbSNP,
					      -dbCore => $dbCore,
					      -dbVariation => $dbVar,
					      -snp_dbname => $sdba->dbc->dbname,
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
					  -snp_dbname => $sdba->dbc->dbname,
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
				      -snp_dbname => $sdba->dbc->dbname,
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
				      -snp_dbname => $sdba->dbc->dbname,
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
				-snp_dbname => $sdba->dbc->dbname,
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

    #and finally, add the alleles to the Allele table for the strains
    #first, insert the single_bp alleles
    $dbVariation->do(qq{INSERT INTO allele (variation_id, sample_id, allele, frequency)
				      SELECT ig.variation_id, ip.population_sample_id, ig.allele_1, 1
				      FROM individual_genotype_single_bp ig, individual_population ip, sample s
				      WHERE ig.sample_id = ip.individual_sample_id
				      AND   ip.population_sample_id = s.sample_id
				      AND   s.is_strain = 1
				  });
    #do the same for the multiple_bp alleles
    $dbVariation->do(qq{INSERT INTO allele (variation_id, sample_id, allele, frequency)
				      SELECT ig.variation_id, ip.population_sample_id, ig.allele_1, 1
				      FROM individual_genotype_multiple_bp ig, individual_population ip, sample s
				      WHERE ig.sample_id = ip.individual_sample_id
				      AND   ip.population_sample_id = s.sample_id
				      AND   s.is_strain = 1
				  });
    #and remove the extra column in sample to identify strains
    $dbVariation->do(qq{ALTER TABLE sample DROP COLUMN is_strain});
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
      -mapping_file <filename>  file containing the mapping data
EOF

      die("\n$msg\n\n");
}

