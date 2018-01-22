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
use Getopt::Long;
use DBI qw(:sql_types);
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::OntologyDBAdaptor;
use ImportUtils;
use String::Approx qw(amatch adist);
use Algorithm::Diff qw(diff);
use XML::LibXML;
use Text::CSV;
use HTTP::Tiny;
use JSON;
use FileHandle;
# set output autoflush for progress bars
$| = 1;

my (
  $infile, $help,
  $host, $dbname, $user, $pass, $port,
  $chost, $cdbname, $cuser, $cpass, $cport,
  $ontology_host, $ontology_dbname, $ontology_user, $ontology_port,
  $skip_synonyms, $skip_phenotypes, $skip_sets,
  $source, $source_version, $set, $assembly, $threshold,
  $global_phenotype_description, $mim2gene, $working_dir, $coord_file,
);

our $verbose;

my %SOURCES = (
  "Uniprot" => {
    description => "Variants with protein annotation imported from Uniprot",
    url => "http://www.uniprot.org/",
    set => "ph_uniprot",
    type => "Variation",
  },
  
  "OMIM" => {
    description => "Variations linked to entries in the Online Mendelian Inheritance in Man (OMIM) database",
    url => "http://www.omim.org/",
    set => "ph_omim",
    type => "Variation",
  },
  
  "NHGRI-EBI GWAS catalog" => {
    description => "Variants associated with phenotype data from the NHGRI-EBI GWAS catalog",
    url => "http://www.ebi.ac.uk/gwas/",
    set => "ph_nhgri",
    type => "Variation",
  },
  
  "EGA" => {
    description => "Variants imported from the European Genome-phenome Archive with phenotype association",
    url => "http://www.ebi.ac.uk/ega/",
    type => "Variation",
  },
  
  "OMIA" => {
    description => "Online Mendelian Inheritance in Animals (OMIA) is a database of genes, inherited disorders and traits in more than 200 animal species (other than human and mouse)",
    url => "http://omia.angis.org.au/home/",
    type => "Gene",
  },
  
  "Animal_QTLdb" => {
    description => "The Animal Quantitative Trait Loci (QTL) database (Animal QTLdb) is designed to house all publicly available QTL and association data on livestock animal species",
    url => "http://www.animalgenome.org/cgi-bin/QTLdb/index",
    type => "QTL",
  },
  
  "RGD" => {
    description => "QTLs from the Rat Genome Database (RGD)",
    url => "http://rgd.mcw.edu/",
    type => "QTL",
  },

  "GIANT" => {
    description => "The Genetic Investigation of ANthropometric Traits (GIANT) consortium is an international collaboration that seeks to identify genetic loci that modulate human body size and shape, including height and measures of obesity",
    url => "http://www.broadinstitute.org/collaboration/giant/index.php/Main_Page",
    type => "Variation",
  },
  
  "MAGIC" => {
    description => "MAGIC (the Meta-Analyses of Glucose and Insulin-related traits Consortium) represents a collaborative effort to combine data from multiple GWAS to identify additional loci that impact on glycemic and metabolic traits",
    url => "http://www.magicinvestigators.org/",
    type => "Variation",
  },
  
  "Orphanet" => {
    description => "The portal for rare diseases and drugs",
    url => "http://www.orpha.net/",
    type => "Gene",
  },

  "DDG2P" => {
    description => "Developmental Disorders Genotype-to-Phenotype Database",
    url => "http://decipher.sanger.ac.uk/",
    type => "Gene",
  },
  
  "MIM morbid" => {
    description => "Online Mendelian Inheritance in Man (OMIM) database",
    url => "http://www.omim.org/",
    type => "Gene",
  },

  "dbGaP" => {
    description => "The database of Genotypes and Phenotypes.",
    url => "http://www.ncbi.nlm.nih.gov/gap",
    type => "Variation",
  },
  
  "ZFIN" => {
    description => "The zebrafish model organism database",
    url => "http://www.zfin.org/",
    type => "Gene",
  },
  
  "GOA" => {
    description => "Inferred from Mutant Phenotype (IMP), annotated by Gene Ontology Annotation",
    url => "http://www.ebi.ac.uk/GOA",
    type => "Gene",
    status => "mixed"
  },

  "Cancer Gene Census" => {
    description => "Catalog of genes of which mutations have been causally implicated in cancer",
    url => "http://cancer.sanger.ac.uk/census",
    type => "Gene",
    status => "somatic"
  },
  "IMPC" => {
    description => "International Mouse Phenotyping Consortium",
    url => "http://www.mousephenotype.org/",
    type => "Gene",  
  },
  "MGP" => {
    description => "WTSI Mouse Genetics Project",
    url => "http://www.sanger.ac.uk/resources/mouse/genomes/",
    type => "Gene",  
  },
  "3i" => {
    description => "Infection, Immunity and Immunophenotyping consortium",
    url => "http://www.immunophenotyping.org/",
    type => "Gene",  
  },
  "MGI" => {
    description => "Mouse Genome Informatics",
    url => "http://www.informatics.jax.org/",
    type => "Gene",  
  },

);

usage() if (!scalar(@ARGV));
 
GetOptions(
    'infile=s' => \$infile,
    'host=s' => \$host,
    'dbname=s' => \$dbname,
    'user=s' => \$user,
    'pass=s' => \$pass,
    'port=s' => \$port,
    'cdbname=s' => \$cdbname,
    'chost=s' => \$chost,
    'cuser=s' => \$cuser,
    'cpass=s' => \$cpass,
    'cport=s' => \$cport,
    'ontology_dbname=s' => \$ontology_dbname,
    'ontology_host=s' => \$ontology_host,
    'ontology_user=s' => \$ontology_user,
    'ontology_port=s' => \$ontology_port,
    'assembly=s' => \$assembly,
    'threshold=s' => \$threshold,
    'phenotype=s' => \$global_phenotype_description,
    'source=s' => \$source,
    'version=i' => \$source_version,
    'verbose!' => \$verbose,
    'skip_synonyms!' => \$skip_synonyms,
    'skip_phenotypes!' => \$skip_phenotypes,
    'skip_sets!' => \$skip_sets,
    'help!' => \$help,
    'mim2gene=s' => \$mim2gene,
    'working_dir=s' => \$working_dir,
    'coord_file=s' => \$coord_file,
) or die "ERROR: Failed to parse command line arguments\n";

usage() if ($help);

die ("Database credentials (--host, --dbname, --user, --pass, --port) are required") unless (defined($host) && defined($dbname) && defined($user) && defined($pass));

if ($source =~ /animal.*qtl.*/i) {
  die ("Database credentials (--ontology_host, --ontology_dbname, --ontology_user, --ontology_port) are required for importing AnimalQTL data") unless (defined($ontology_host) && defined($ontology_dbname) && defined($ontology_user));
}

if (! ($source =~ m/^impc/i || ($source =~ m/^mgi/i) )) {
  die ("An input file (--infile) is required") unless (defined($infile));
  die ("A source (--source) is required") unless (defined($source));
  die ("A source version (--version) is required") unless (defined($source_version));
}
$port  ||= 3306;
$cport ||= $port;
$chost ||= $host;
$cuser ||= $user;
#$cpass ||= $pass;
$cpass = '';

my $result;
my $source_name;
my $source_description;
my $source_url;
my $source_status = 'germline';
my $object_type;
my $prev_prog;

my $pubmed_prefix = 'PMID:';

my %special_characters = (
  'ö' => 'o',
  'ü' => 'u',
  'ä' => 'a',
  'í' => 'i',
);



=head

    The parser subroutines parse the input file into a common data structure.
    Currently what is returned should be a reference to a hash. The hash should
    contain the key 'phenotypes' with the value being a reference to an array of
    phenotype data objects. The phenotype data object is a reference to a hash
    where the keys correspond to the column names in the phenotype_feature table.
    In addition, there are the keys 'rsid' which holds the rs-id that the phenotype
    annotates and 'description' and 'name' which correspond to the columns in the
    phenotype table.
    
    In addition, the hash can also contain the key 'synonyms' and the value is a
    reference to a hash where the keys are rs-ids and each respective value is a
    reference to an array of synonyms for the rs-id.

=cut

# Make sure that the input file is XML compliant
#ImportUtils::make_xml_compliant($infile) unless $infile =~ /(gz|z)$/i;

# Remove carriage return in the input file
#remove_carriage_return($infile) unless $infile =~ /(gz|z)$/i;


# Connect to the variation database
print STDOUT localtime() . "\tConnecting to database $dbname\n" if ($verbose);
my $db_adaptor = new Bio::EnsEMBL::Variation::DBSQL::DBAdaptor(
  -host => $host,
  -user => $user,
  -pass => $pass,
  -port => $port,
  -dbname => $dbname
) or die("Could not get a database adaptor for $dbname on $host:$port");
print STDOUT localtime() . "\tConnected to $dbname on $host:$port\n" if ($verbose);

# connect to core DB for genes
my $core_db_adaptor;

if(defined($chost) && defined($cuser) && defined($cdbname)) {
  print STDOUT localtime() . "\tConnecting to database $cdbname\n" if ($verbose);
  $core_db_adaptor = new Bio::EnsEMBL::DBSQL::DBAdaptor(
    -host => $chost,
    -user => $cuser,
    -pass => $cpass,
    -port => $cport,
    -dbname => $cdbname
  ) or die("Could not get a database adaptor for $cdbname on $chost:$cport");
  print STDOUT localtime() . "\tConnected to $cdbname on $chost:$cport\n" if ($verbose);
}

# connect to ontology DB for mapping CMO names to accessions for AnimalQTL import 
my $ontology_db_adaptor;

if(defined($ontology_host) && defined($ontology_user) && defined($ontology_dbname)) {
  print STDOUT localtime() . "\tConnecting to database $ontology_dbname\n" if ($verbose);
  $ontology_db_adaptor = new Bio::EnsEMBL::DBSQL::OntologyDBAdaptor(
    -host => $ontology_host,
    -user => $ontology_user,
    -port => $ontology_port,
    -dbname => $ontology_dbname
  ) or die("Could not get a database adaptor for $ontology_dbname on $ontology_host:$ontology_port");
  print STDOUT localtime() . "\tConnected to $ontology_dbname on $ontology_host:$ontology_port\n" if ($verbose);
}

# get attrib types
our @attrib_types = @{get_attrib_types($db_adaptor)};

# get phenotypes from DB
my $phenotype_adaptor = $db_adaptor->get_PhenotypeAdaptor;
our %phenotype_cache = map {$_->description() => $_->dbID()} @{$phenotype_adaptor->fetch_all};

my $initial_phenotype_count = $phenotype_adaptor->generic_count;

print STDOUT localtime() . "\tParsing input file\n" if ($verbose);

my @ids;

# Parse the input files into a hash
if ($source =~ m/uniprot/i) {
  $skip_phenotypes = 1;
  $result = parse_uniprot($infile);
  $source_name = 'Uniprot';
}
elsif ($source =~ m/nhgri/i) {
  $result = parse_nhgri($infile);
  $source_name = 'NHGRI-EBI GWAS catalog';
}
elsif($source =~ /omim.*gene/i) {
  #die("ERROR: No core DB parameters supplied (--chost, --cdbname, --cuser) or could not connect to core database") unless defined($core_db_adaptor);
  #die("ERROR: mim2gene file required (--mim2gene [file])\n") unless defined($mim2gene) && -e $mim2gene;
  
  #$result = parse_omim_gene($infile, $mim2gene, $core_db_adaptor);
  $result = parse_omim_gene($infile);
  $source_name = 'MIM morbid';
}
elsif ($source =~ m/^omim$/i) {
  $result = parse_omim($infile);
  $source_name = 'OMIM';
}
elsif ($source =~ m/ega/i) {
  $source_name = 'EGA';
  $source_description = $SOURCES{$source_name}->{description};
  $source_url         = $SOURCES{$source_name}->{url};
  my $source_id = get_or_add_source($source_name,$source_description,$source_url,$source_status,$db_adaptor);
  print STDOUT "$source source_id is $source_id\n" if ($verbose);
  parse_ega($infile,$source_id);
  exit(0);
}
elsif ($source =~ /omia/i) {
  die("ERROR: No core DB parameters supplied (--chost, --cdbname, --cuser) or could not connect to core database") unless defined($core_db_adaptor); 
  
  $result = parse_omia($infile, $core_db_adaptor);
  $source_name = 'OMIA';
}
elsif($source =~ /animal.*qtl.*/i) {
  $result = parse_animal_qtl($infile, $db_adaptor);
  $source_name = 'Animal_QTLdb';
}
elsif($source =~ /rgd_qtl/i) {
  die("ERROR: Assembly version not specified - use --assembly [version]") unless defined($assembly);
  
  $result = parse_rgd_qtl($infile, $assembly, $db_adaptor);
  $source_name = 'RGD';
}
elsif($source =~ /rgd_gene/i) {
  die("ERROR: No core DB parameters supplied (--chost, --cdbname, --cuser) or could not connect to core database") unless defined($core_db_adaptor);

  $result = parse_rgd_gene($infile, $db_adaptor);
  $source_name = 'RGD';
  $SOURCES{$source_name}->{type} = 'Gene'; # By default it is set to 'QTL'
}
elsif($source =~ /giant/i) {
  die("ERROR: p-value threshold not specified - use --threshold [p-value]") unless defined($threshold);
  die("ERROR: phenotype description not specified - use --phenotype [description]") unless defined($global_phenotype_description);
  
  $result = parse_giant($infile);
  $source_name = 'GIANT';
}
elsif($source =~ /magic/i) {
  die("ERROR: p-value threshold not specified - use --threshold [p-value]") unless defined($threshold);
  die("ERROR: phenotype description not specified - use --phenotype [description]") unless defined($global_phenotype_description);
  
  $result = parse_magic($infile);
  $source_name = 'MAGIC';
}
elsif($source =~ /orphanet/i) {
  die("ERROR: No core DB parameters supplied (--chost, --cdbname, --cuser) or could not connect to core database") unless defined($core_db_adaptor);
  
  $result = parse_orphanet($infile, $core_db_adaptor);
  $source_name = 'Orphanet';
}
elsif($source =~ /ddg2p/i) {
  die("ERROR: No core DB parameters supplied (--chost, --cdbname, --cuser) or could not connect to core database") unless defined($core_db_adaptor);
  
  $result = parse_ddg2p($infile, $core_db_adaptor);
  $source_name = 'DDG2P';
}
elsif($source =~ /^cgc$/i) {
  die("ERROR: No core DB parameters supplied (--chost, --cdbname, --cuser) or could not connect to core database") unless defined($core_db_adaptor);

  $result = parse_cancer_gene_census($infile, $core_db_adaptor);
  $source_name = 'Cancer Gene Census';
}
elsif($source =~ /mim.+dump/) {
  $result = parse_mim_dump($infile);
  $source_name = 'MIM morbid';
}
elsif ($source =~ m/zfin/i) {
  $result = parse_zfin($infile, $core_db_adaptor);
  $source_name = 'ZFIN';
}
elsif ($source =~ m/dbgap/i) {
  $result = parse_dbgap($infile);
  $source_name = 'dbGaP';
}
elsif ($source =~ m/^goa$/i) {
  $result = parse_goa($infile, $core_db_adaptor);
  $source_name = 'GOA';
}
elsif ($source =~ m/^impc/i || $source =~ m/^mgi/i) {
  $skip_synonyms = 1;
  my ($data_source, $url);
  if ($source =~ m/^impc/i ) {
    $data_source = 'impc';
    $source_name = 'IMPC';
    $url = '/mi/impc/solr/genotype-phenotype';
  } else {
    $data_source = 'mgi';
    $source_name = 'MGI';
    $url = '/mi/impc/solr/mgi-phenotype';
  }
  die "Coord file is required: ftp://ftp.informatics.jax.org/pub/reports/MGI_MRK_Coord.rpt" if (!-f $coord_file);
  if (!$infile) {
    $infile = get_mouse_phenotype_data($working_dir, $data_source, $url);
  }
  my $mouse_phenotype_source_ids = get_mouse_phenotype_data_source_ids($infile, $data_source);
  $source_version = update_mouse_phenotype_data_version($mouse_phenotype_source_ids); 
  clear_mouse_phenotype_data_from_last_release($mouse_phenotype_source_ids);
  my $marker_coords = get_marker_coords($infile, $coord_file);
  $result = parse_mouse_phenotype_data($infile, $marker_coords, $data_source, $mouse_phenotype_source_ids);
}

else {
  die("Source $source is not recognized");
}

$source_description = $SOURCES{$source_name}->{description};
$source_url         = $SOURCES{$source_name}->{url};
$set                = defined($SOURCES{$source_name}->{set}) ? $SOURCES{$source_name}->{set} : undef;
$object_type        = $SOURCES{$source_name}->{type};
$source_status      = $SOURCES{$source_name}->{status} if (defined($SOURCES{$source_name}->{status}));

my %synonym;
my @phenotypes;
if (exists($result->{'synonyms'})) {
  %synonym = %{$result->{'synonyms'}};
  # To get all the ids of the source (Uniprot)
  @ids = keys(%synonym);
}
if (exists($result->{'phenotypes'})) {
  @phenotypes = @{$result->{'phenotypes'}};
}

print STDOUT "Got ".(scalar @phenotypes)." objects\n" if ($verbose);

# Get internal variation ids for the rsIds
print STDOUT "Retrieving internal variation IDs\n" if ($verbose);
if (scalar @ids == 0) {
  @ids = map {$_->{'id'}} @phenotypes;
}
my $variation_ids = $object_type =~ /Variation/ ? get_dbIDs(\@ids,$db_adaptor) : {};

# Get coordinates of objects
my $coords;

print STDOUT "Retrieving object coordinates\n" if ($verbose);

# might be able to copy them from data (QTLs and Genes)
if(defined($phenotypes[0]->{seq_region_id})) {
  foreach my $p(@phenotypes) {
    my $coord = {};
    
    foreach my $key(qw(id start end strand)) {
      $coord->{'seq_region_'.$key} = $p->{'seq_region_'.$key};
    }
    
    push @{$coords->{$p->{id}}}, $coord;
  }
}

$coords ||= get_coords(\@ids, $object_type, $object_type eq 'Gene' ? $core_db_adaptor : $db_adaptor);
  
# uniquify coords
foreach my $id (keys %$coords) {
  my %tmp = map {
    $_->{seq_region_id}."_".
    $_->{seq_region_start}."_".
    $_->{seq_region_end}."_".
    $_->{seq_region_strand}."_" => $_
  } @{$coords->{$id}};
  
  $coords->{$id} = [values %tmp];
}

# Get or add a source
my $source_id = get_or_add_source($source_name,$source_description,$source_url,$source_status,$db_adaptor);
print STDOUT "$source source_id is $source_id\n" if ($verbose);


# Add the synonyms if required
unless ($skip_synonyms) {
  print STDOUT "Adding synonyms\n" if ($verbose);
  add_synonyms(\%synonym,$variation_ids,$source_id,$db_adaptor);
}

# Now, insert phenotypes
unless ($skip_phenotypes) {
  die("ERROR: No phenotypes or objects retrieved from input\n") unless scalar @phenotypes;

  print STDOUT "Adding phenotypes\n" if ($verbose);
  add_phenotypes(\@phenotypes,$coords,$source_id,$object_type,$db_adaptor);

  my $added_phenotypes = $phenotype_adaptor->generic_count - $initial_phenotype_count;
  print STDOUT "$added_phenotypes new phenotypes added\n" if ($verbose);
}

# Add the variation sets if required
unless ($skip_sets) {
  if (%synonym && !$skip_synonyms && $skip_phenotypes) {
    print STDOUT "Adding variation sets for synonyms\n" if ($verbose);
    add_set($set,$source_id,$db_adaptor,'synonym');
  }
  elsif (@phenotypes && !$skip_phenotypes) {
    print STDOUT "Adding variation sets for phenotypes\n" if ($verbose);
    add_set($set,$source_id,$db_adaptor,'phenotype');
  }
}



# Loop over the remaining ids (the ones that could not be find in the db) and print them out
#while (my ($rs_id,$var_id) = each(%{$variation_ids})) {
#    next if (defined($var_id->[0]));
#    print STDOUT "$rs_id could not be found in $dbname";
#    if (defined($synonym{$rs_id})) {
#        print STDOUT " (Synonyms: " . join(", ",@{$synonym{$rs_id}}) . ")";
#    }
#    print STDOUT "\n";
#}



###########
# METHODS #
###########

sub parse_uniprot {
  my $infile = shift;
  
  my %synonym;
  my @phenotypes;
  
  # Open the input file for reading
  if($infile =~ /gz$/) {
    open IN, "zcat $infile |" or die ("Could not open $infile for reading");
  }
  else {
    open(IN,'<',$infile) or die ("Could not open $infile for reading");
  }
  
  # Read through the file and parse out the desired fields
  while (<IN>) {
    chomp;
    
    # A regexp to catch the meta information in the header. Just echo this to the stdout for logging purposes
    if ($_ =~ m/^\s*(Description|Name|Release)\:/) {
      print STDOUT $_ . "\n";
    }
    
    # Main regexp to extract relevant variation information
    if ($_ =~ m/^(\S+)\s+\S+\s+(VAR\_\d+)\s+\w\.\S+\s+(Disease|Polymorphism|Unclassified)\s+(\-|rs\d*)\s*(.*)$/) {
      
      # Get the data that was caught by the regexp
      my $gene = $1;
      my $uniprot_id = $2;
      my $rs_id = $4;
      my $phenotype = $5;
      
      # If no rsId was given, will attempt to get one by looking up the Uniprot id in the synonym table
      if ($rs_id ne '-') {
        push(@{$synonym{$rs_id}},$uniprot_id);
      }
      #else {
      #  $rs_id = $uniprot_id;
      #}
      
      #$phenotype ||= '-';
      
      ## Try to further split the phenotype into a short name, description and possibly MIM id
      #if ($phenotype ne '-') {
      #  my $description;
      #  my $name;
      #  my $mim_id;
      #  
      #  ($description,$mim_id) = $phenotype =~ m/^([^\[]+)(?:\[(MIM\:.+?)\])?$/;
      #  ($description,$name) = $description =~ m/^(.+?)\s*(?:\((.+?)\))?\s*$/;
      #  
      #  $mim_id &&= join(",MIM:",split(",",$mim_id));
      #  $mim_id =~ s/\s//g if (defined($mim_id));
      #
      #  push(
      #    @phenotypes,
      #    {
      #      "id"              => $rs_id,
      #      "associated_gene" => $gene,
      #      "description"     => $description,
      #      "name"            => $name,
      #      "variation_names" => $rs_id,
      #      "study"           => $mim_id,
      #      "study_type"      => 'GWAS'
      #    }
      #  );
      #}
    }
  }
  close(IN);
  
  print STDOUT "Parsed " . scalar(keys(%synonym)) . " rs-ids with Uniprot synonyms\n" if ($verbose);
  #print STDOUT scalar(@phenotypes) . " phenotype associations were found linked to rs-ids\n" if ($verbose);
  
  #my %result = ('synonyms' => \%synonym, 'phenotypes' => \@phenotypes);
  my %result = ('synonyms' => \%synonym);
  return \%result;
}

sub parse_nhgri {
  my $infile = shift;
  
  my %headers;
  my @phenotypes;
  
  # Open the input file for reading
  if($infile =~ /gz$/) {
    open IN, "zcat $infile |" or die ("Could not open $infile for reading");
  }
  else {
    open(IN,'<',$infile) or die ("Could not open $infile for reading");
  }
  
  # Read through the file and parse out the desired fields
  while (<IN>) {
    chomp;

    my @row_data = split(/\t/,$_);
    
    # header
    if(/^DATE\s+ADDED\s+TO\s+CATALOG/) {
      $headers{uc($row_data[$_])} = $_ for 0..$#row_data;
    }
    else {
      die "ERROR: Couldn't find header data\n" unless %headers;
    
      my %content;
      $content{$_} = $row_data[$headers{$_}] for keys %headers;
    
      my $pubmed_id      = $content{'PUBMEDID'};
      my $study          = $content{'STUDY'};
      my $phenotype      = $content{'DISEASE/TRAIT'};
      my $gene           = ($content{'REPORTED GENE(S)'} =~ /\?/) ? '' : $content{'REPORTED GENE(S)'};
      my $rs_risk_allele = ($content{'STRONGEST SNP-RISK ALLELE'} =~ /\?/) ? '' : $content{'STRONGEST SNP-RISK ALLELE'};
      my $rs_id          = $content{'SNPS'};
      my $pvalue         = ($content{'P-VALUE'} ne '') ? $content{'P-VALUE'} : '';
      my $ratio          = $content{'OR OR BETA'};
      my $ratio_info     = $content{'95% CI (TEXT)'};
      my @accessions     = split/\,/, $content{'MAPPED_TRAIT_URI'};

      my $risk_frequency = '';
      if ($rs_risk_allele =~ /^\s*$rs_id-+\s*(\w+)\s*$/i) {
        $rs_risk_allele = $1;
        $risk_frequency = $content{'RISK ALLELE FREQUENCY'};
      }
      
      $gene =~ s/\s+//g;
      $gene =~ s/–/-/g;

      my %data = (
        'study_type' => 'GWAS',
        'description' => $phenotype,
        'associated_gene' => $gene,
        'risk_allele' => $rs_risk_allele,
        'risk_allele_freq_in_controls' => $risk_frequency,
        'p_value' => $pvalue,
        'study_description' => $study,
        'accessions'   => \@accessions,
        'ontology_mapping_type' =>'is'
      );
   
    
      # Post process the ratio data
      if (defined($ratio)) {
        if ($ratio =~ /(\d+)?(\.\d+)$/) {
          my $pre  = $1;
          my $post = $2;
          $ratio = (defined($pre)) ? "$pre$post" : "0$post";
          $ratio = 0 if ($ratio eq '0.00');
        } else {
          $ratio = undef;
        }
      }

      # Add ratio/coef
      if (defined($ratio)) {
        # Parse the ratio info column to extract the unit information (we are not interested in the confidence interval)
        if ($ratio_info =~ /^\s*(\[.+\])?\s*(.+)$/) {
          my $unit = $2;
             $unit =~ s/\(//g;
             $unit =~ s/\)//g;
             $unit =~ s/µ/micro/g;
          if ($unit =~ /^\s+$/ || $ratio >= 1) {
            $data{'odds_ratio'} = $ratio;
          }
          else {
            $data{'beta_coef'} = "$ratio $unit";
          }
        }
        else {
          $data{'odds_ratio'} = $ratio;
        }
      }
    
      # Parse the ids
      my @ids;
      $rs_id ||= "";
      while ($rs_id =~ m/(rs[0-9]+)/g) {
        push(@ids,$1);
      }
      $data{'variation_names'} = join(',',@ids);
      $data{'study'} = $pubmed_prefix . $pubmed_id if (defined($pubmed_id));
      
      # If we didn't get any rsIds, skip this row (this will also get rid of the header)
      warn("Could not parse any rsIds from string '$rs_id'") if (!scalar(@ids));
      next if (!scalar(@ids));
      
      map {
        my %t_data = %{\%data};
        $t_data{'id'} = $_;
        push(@phenotypes,\%t_data)
      } @ids;
    }
  }
  close(IN);

  my %result = ('phenotypes' => \@phenotypes);
  
  return \%result;
}

sub parse_omim {
  my $infile = shift;
  
  my @phenotypes;
  
  # Open the input file for reading
  if($infile =~ /gz$/) {
    open IN, "zcat $infile |" or die ("Could not open $infile for reading");
  }
  else {
    open(IN,'<',$infile) or die ("Could not open $infile for reading");
  }
  
  # Read through the file and parse out the desired fields
  while (<IN>) {
    chomp;
    
    my @attributes = split(/\t/);
    
    next if (!$attributes[4] or $attributes[0] =~ /^\D/);
    
    my ($study,$allele) = split(/\./,$attributes[0]);
    my $ids = $attributes[4];
    
    my @rs_list = split(',',$ids);
    
    #  Get one line for each variation_names
    foreach my $rs (@rs_list) {

      # Skip the risk allele if the variant is "0000"
      my $data = {
        'id'      => $rs,
        'study'       => 'MIM:'.$study,
        'associated_gene' => $attributes[2],
        'variation_names' => $ids,
        'description'   => $attributes[1],
        'risk_allele' => ($allele !~ m/^\s*0+\s*$/ ? $allele : undef),
      };
    
      push(@phenotypes,$data);
    }
  }
  close(IN);
  
  my %result = ('phenotypes' => \@phenotypes);
  return \%result;
}

sub parse_ega {
  my $infile = shift;
  my $source_id = shift;
  
  my $study_check_stmt = qq{
    SELECT
      study_id
    FROM
      study
    WHERE
      name=? AND source_id=$source_id
    LIMIT 1
  };
  my $nhgri_check_stmt = qq{
    SELECT
      st.study_id,st.study_type
    FROM
      study st, source s
    WHERE
      external_reference=? AND s.name like '%nhgri%'
      AND s.source_id=st.source_id
    LIMIT 1
  };
  my $study_ins_stmt = qq{
    INSERT INTO
      study (
      name,
      source_id,
      external_reference,
      url,
      study_type
      )
    VALUES (
      ?,
      $source_id,
      ?,
      ?,
      ?
    )
  };
  # NHGRI and EGA associated studies
  my $asso_study_check_stmt = qq{
    SELECT
      study1_id
    FROM
      associate_study
    WHERE
      study1_id = ? AND study2_id = ?
    LIMIT 1
  };
  my $asso_study_ins_stmt = qq{
    INSERT INTO
      associate_study (study1_id,study2_id)
    VALUES (?,?)
  };
  
  my $nhgri_check_sth = $db_adaptor->dbc->prepare($nhgri_check_stmt);
  my $study_check_sth = $db_adaptor->dbc->prepare($study_check_stmt);
  my $study_ins_sth   = $db_adaptor->dbc->prepare($study_ins_stmt);
  my $asso_study_check_sth = $db_adaptor->dbc->prepare($asso_study_check_stmt);
  my $asso_study_ins_sth   = $db_adaptor->dbc->prepare($asso_study_ins_stmt);
  
  # Open the input file for reading
  open(IN,'<',$infile) or die ("Could not open $infile for reading");
  
  # Read through the file and parse out the desired fields
  while (<IN>) {
    chomp $_;
    my @attributes = split(",",$_);
    next if ($attributes[1] eq '');
    my $name = $attributes[0];
    my $pubmed = $pubmed_prefix.$attributes[1];
    my $url = $attributes[2];
    
    # NHGRI study
    my $nhgri_study_id;
    my $study_type;
    $nhgri_check_sth->bind_param(1,$pubmed,SQL_VARCHAR);
    $nhgri_check_sth->execute();
    $nhgri_check_sth->bind_columns(\$nhgri_study_id,\$study_type);
    $nhgri_check_sth->fetch();
    
    if (!defined($nhgri_study_id)) {
      print "No NHGRI-EBI study found for the EGA $name | $pubmed !\n";
      next;
    }
    
    # EGA study
    my $study_id;
    $study_check_sth->bind_param(1,$name,SQL_VARCHAR);
    $study_check_sth->execute();
    $study_check_sth->bind_columns(\$study_id);
    $study_check_sth->fetch();
    if (!defined($study_id)) {
      $study_ins_sth->bind_param(1,$name,SQL_VARCHAR);
      $study_ins_sth->bind_param(2,$pubmed,SQL_VARCHAR);
      $study_ins_sth->bind_param(3,$url,SQL_VARCHAR);
      $study_ins_sth->bind_param(4,$study_type,SQL_VARCHAR);
      $study_ins_sth->execute();
      $study_id = $db_adaptor->dbc->db_handle->{'mysql_insertid'};
    }
    
    my $is_associated;
    $asso_study_check_sth->bind_param(1,$nhgri_study_id,SQL_INTEGER);
    $asso_study_check_sth->bind_param(2,$study_id,SQL_INTEGER);
    $asso_study_check_sth->execute();
    $asso_study_check_sth->bind_columns(\$is_associated);
    $asso_study_check_sth->fetch();
    
    if (!defined($is_associated)) {
      $asso_study_ins_sth->bind_param(1,$nhgri_study_id,SQL_INTEGER);
      $asso_study_ins_sth->bind_param(2,$study_id,SQL_INTEGER);
      $asso_study_ins_sth->execute();
    }
  }
  close(IN);
}

sub parse_goa {
  my $infile = shift;
  my $core_db_adaptor = shift;
  
  my $ga = $core_db_adaptor->get_GeneAdaptor;
  die("ERROR: Could not get gene adaptor") unless defined($ga);
  
  my @phenotypes;
  my %headers;
  
  # Open the input file for reading
  if($infile =~ /gz$/) {
    open IN, "zcat $infile |" or die ("Could not open $infile for reading");
  }
  else {
    open(IN,'<',$infile) or die ("Could not open $infile for reading");
  }
  
  # Read through the file and parse out the desired fields
  while (<IN>) {
    chomp;
    
    my @row_data = split(/\t/,$_);
    
    # header
    if(/^Gene_ID/) {
      $headers{uc($row_data[$_])} = $_ for 0..$#row_data;
    }
    else {
      die "ERROR: Couldn't find header data\n" unless %headers;
    
      my %content;
      $content{$_} = $row_data[$headers{$_}] for keys %headers;

      my $pubmed_ids  = ($content{'PMIDS'} eq '') ? undef : $content{'PMIDS'};
      my $gene_id     = $content{'GENE_ID'};
      my $external_id = $content{'GO_TERM'};
      my $phenotype   = $content{'GO_DESCRIPTION'};
      my $uniprot_ids = $content{'UNIPROT_IDS'};
     
    
      my $gene = $ga->fetch_by_stable_id($gene_id);
    
      if(!$gene) {
        print STDERR "WARNING: No matching found the Ensembl gene ID '$gene_id'\n";
      }
    
      next unless $gene;
    
    
      my %data = (
        'id' => $gene->stable_id,
        'description' => $phenotype,
        'external_id' => $external_id,
        'seq_region_id' => $gene->slice->get_seq_region_id,
        'seq_region_start' => $gene->seq_region_start,
        'seq_region_end' => $gene->seq_region_end,
        'seq_region_strand' => $gene->seq_region_strand,
        'xref_id' => $uniprot_ids,
        'accessions' => [$external_id],
         ontology_mapping_type =>'is' 
      );
      
      $data{'study'} = $pubmed_ids if ($pubmed_ids);
    
      push(@phenotypes,\%data);
    }
  }
  close(IN);
  
  my %result = ('phenotypes' => \@phenotypes);
  return \%result;
}

sub parse_omia {
  my $infile = shift;
  my $core_db_adaptor = shift;
  
  my $ga = $core_db_adaptor->get_GeneAdaptor;
  die("ERROR: Could not get gene adaptor") unless defined($ga);
  
  my @phenotypes;
  
  # Open the input file for reading
  if($infile =~ /gz$/) {
    open IN, "zcat $infile |" or die ("Could not open $infile for reading");
  }
  else {
    open(IN,'<',$infile) or die ("Could not open $infile for reading");
  }
  
  # Read through the file and parse out the desired fields
  while (<IN>) {
    chomp;
    
    next if /^gene_symbol/;
    
    my @data = split /\t/, $_;
    
    my $genes = $ga->fetch_all_by_external_name($data[0]);
    
    if(scalar @$genes != 1) {
      print STDERR "WARNING: Found ".(scalar @$genes)." matching Ensembl genes for gene name $data[0]\n";
    }
    
    next unless scalar @$genes;
    
    foreach my $gene(@$genes) {
      push @phenotypes, {
        'id' => $gene->stable_id,
        'description' => $data[5],
        'external_id' => 'OMIA'.$data[2],
        'seq_region_id' => $gene->slice->get_seq_region_id,
        'seq_region_start' => $gene->seq_region_start,
        'seq_region_end' => $gene->seq_region_end,
        'seq_region_strand' => $gene->seq_region_strand
      };
    }
  }
  close(IN);
  
  my %result = ('phenotypes' => \@phenotypes);
  return \%result;
}

sub parse_animal_qtl {
  my $infile = shift;
  my $db_adaptor = shift;
 
  my $ontology_term_adaptor = $ontology_db_adaptor->get_OntologyTermAdaptor; 
  # get seq_region_ids
  my $seq_region_ids = get_seq_region_ids($db_adaptor);
  
  my @phenotypes;
  
  # Open the input file for reading
  if($infile =~ /gz$/) {
    open IN, "zcat $infile |" or die ("Could not open $infile for reading");
  }
  else {
    open(IN,'<',$infile) or die ("Could not open $infile for reading");
  }
  
  # Read through the file and parse out the desired fields
  while (<IN>) {
    chomp;
    
    next if /^(\#|\s)/ || !$_;
    
    my @data = split /\t/, $_;
    
    # fix chr
    $data[0] =~ s/^chr(om)?\.?//i;
    
    if(!defined($seq_region_ids->{$data[0]})) {
      print STDERR "WARNING: Could not find seq_region_id for chromosome name $data[0]\n";
      next;
    }

    # parse "extra" GFF fields
    my $extra = {};
    
    foreach my $chunk(split /\;/, $data[8]) {
      my ($key, $value) = split /\=/, $chunk;
      next unless defined($key) && defined($value);
      $value =~ s/\"|\'//g;
      $extra->{$key} = $value;
    }
    
    if ($extra->{'Map_Type'} eq 'Linkage' || $data[3] == 0) {
      print STDERR "WARNING: Could not find a precise location for the QTL ".$extra->{QTL_ID}."\n";
      next;
    }

    if ($data[4] !~ /^\d+$/) {
      print STDERR "WARNING: Could not find a numeric seq_region_end for the QTL ".$extra->{QTL_ID}."\n";
      next;
    }

    # create phenotype hash

    my $phenotype = {
      'id' => $extra->{QTL_ID},
      'description' => $extra->{Name},
      'seq_region_id' => $seq_region_ids->{$data[0]},
      'seq_region_start' => $data[3] || 1,
      'seq_region_end' => $data[4],
      'seq_region_strand' => 1
    };
    my $CMO_name = $extra->{CMO_name};
    
    if ($CMO_name) {
      my @accessions = ();
      my $terms = $ontology_term_adaptor->fetch_all_by_name($CMO_name, 'CMO'); 
      foreach my $term (@$terms) {
        push @accessions, $term->accession;
      }
      $phenotype->{ontology_mapping_type} = 'is';
      $phenotype->{accessions} = \@accessions;
    }
    
    if ($phenotype->{'seq_region_start'} > $phenotype->{'seq_region_end'}) {
      my $tmp_end = $phenotype->{'seq_region_end'};
      $phenotype->{'seq_region_start'} = $phenotype->{'seq_region_end'};
      $phenotype->{'seq_region_end'}   = $tmp_end;
    }
    
    # add additional fields if found
    $phenotype->{'study'} = $pubmed_prefix.$extra->{'PUBMED_ID'} if defined($extra->{'PUBMED_ID'} && $extra->{'PUBMED_ID'} =~ /^\d+$/);
    $phenotype->{'p_value'} = $extra->{'P-value'} if defined($extra->{'P-value'});
    $phenotype->{'f_stat'} = $extra->{'F-stat'} if defined($extra->{'F-stat'});
    $phenotype->{'lod_score'} = $extra->{'LOD-score'} if defined($extra->{'LOD-score'});
    $phenotype->{'variance'} = $extra->{'Variance'} if defined($extra->{'Variance'});
    $phenotype->{'associated_gene'} = $extra->{'Candidate_Gene_Symble'} if defined($extra->{'Candidate_Gene_Symble'});
    
    push @phenotypes, $phenotype;
  }
  
  my %result = ('phenotypes' => \@phenotypes);
  return \%result;
}

sub parse_rgd_qtl {
  my $infile = shift;
  my $assembly = shift;
  my $db_adaptor = shift;
  
  # get seq_region_ids
  my $seq_region_ids = get_seq_region_ids($db_adaptor);
  
  my @phenotypes;
  
  # Open the input file for reading
  if($infile =~ /gz$/) {
    open IN, "zcat $infile |" or die ("Could not open $infile for reading");
  }
  else {
    open(IN,'<',$infile) or die ("Could not open $infile for reading");
  }
  
  my (%headers, $line_num);
  
  # Read through the file and parse out the desired fields
  while (<IN>) {
    chomp;
    $line_num++;
    
    next if /^(\#|\s)/ || !$_;
    
    my @data = split /\t/, $_;
    
    # header
    if(/^QTL_RGD_ID/) {
      $headers{$data[$_]} = $_ for 0..$#data;
    }
    else {
      die "ERROR: Couldn't find header data\n" unless %headers;
      
      my %data;
      $data{$_} = $data[$headers{$_}] for keys %headers;
      
      # check chromosome
      my $chr = $data{$assembly.'_MAP_POS_CHR'};
      
      if(!defined($chr) || !$chr) {
        print STDERR "WARNING: Could not get coordinates for assembly $assembly on line $line_num\n";
        next;
      }
      
      $chr =~ s/^chr(om)?\.?//i;
      
      if(!defined($seq_region_ids->{$chr})) {
        print STDERR "WARNING: Could not find seq_region_id for chromosome name $chr on line $line_num\n";
        next;
      }
      
      my $description = $data{TRAIT_NAME};
      if ($description =~ /\(/) {
        $description = (split(/\s\(/, $description))[0];
      }
      
      my $phenotype = {
        id => $data{QTL_SYMBOL},
        description => $description,
        seq_region_id => $seq_region_ids->{$chr},
        seq_region_start => $data{$assembly.'_MAP_POS_START'},
        seq_region_end => $data{$assembly.'_MAP_POS_STOP'},
        seq_region_strand => 1,
        lod_score => $data{LOD},
        p_value => $data{P_VALUE},
        variance => $data{VARIANCE},
        associated_gene => $data{CANDIDATE_GENE_SYMBOLS}
      };
      
      push @phenotypes, $phenotype;
    }
  }
  
  my %result = ('phenotypes' => \@phenotypes);
  return \%result;
}

sub parse_rgd_gene {
  my $infile = shift;
  my $db_adaptor = shift;

  # get seq_region_ids
  my $seq_region_ids = get_seq_region_ids($db_adaptor);

  my @phenotypes;

  my $gene_adaptor = $core_db_adaptor->get_GeneAdaptor;
  die("ERROR: Could not get gene adaptor") unless defined($gene_adaptor);

  # Open the input file for reading
  if($infile =~ /gz$/) {
    open IN, "zcat $infile |" or die ("Could not open $infile for reading");
  }
  else {
    open(IN,'<',$infile) or die ("Could not open $infile for reading");
  }

  my %rgd_coords;

  my (%headers, $line_num);

  # Read through the file and parse out the desired fields
  while (<IN>) {
    chomp;
    $line_num++;

    next if /^(\#|\s)/ || !$_;

    my @data_line = split /\t/, $_;

    # header
    if(/^RGD_ID/) {
      $headers{$data_line[$_]} = $_ for 0..$#data_line;
    }
    else {
      die "ERROR: Couldn't find header data\n" unless %headers;

      my %data;
      $data{$_} = $data_line[$headers{$_}] for keys %headers;

      next unless ($data{'OBJECT_TYPE'} =~ /gene/i);

      my $symbol = $data{'OBJECT_SYMBOL'};

      if (!$rgd_coords{$symbol}) {

        my $gene = $gene_adaptor->fetch_by_display_label($symbol);

        if (!$gene) {
          my $genes = $gene_adaptor->fetch_all_by_external_name($symbol,'RGD');
          $gene = $genes->[0] if (scalar(@$genes) > 0);
        }

        if (!$gene) {
          print STDERR "Symbol '$symbol' not found in the Core database\n";
          next;
        }

        $rgd_coords{$symbol} = { 'gene'   => $gene->stable_id,
                                 'chr'    => $gene->slice->seq_region_name,
                                 'start'  => $gene->start,
                                 'end'    => $gene->end,
                                 'strand' => $gene->strand
                               };
      }

      my $phenotype = {
        id => $rgd_coords{$symbol}{'gene'},
        description => $data{'TERM_NAME'},
        seq_region_id => $seq_region_ids->{$rgd_coords{$symbol}{'chr'}},
        seq_region_start => $rgd_coords{$symbol}{'start'},
        seq_region_end => $rgd_coords{$symbol}{'end'},
        seq_region_strand => $rgd_coords{$symbol}{'strand'},
      };

      if ($data{'REFERENCES'} =~ /^(RGD\:\d+)\|/) {
        $phenotype->{external_id} = $1;
      }

      if ($data{'TERM_ACC_ID'} =~ /^MP\:\d+/) {
        $phenotype->{accessions} = [ $data{'TERM_ACC_ID'} ];
      }

      push @phenotypes, $phenotype;
    }
  }

  my %result = ('phenotypes' => \@phenotypes);
  return \%result;
}


sub parse_giant {
  my $infile = shift;
  
  my @phenotypes;
  
  # Open the input file for reading
  if($infile =~ /gz$/) {
    open IN, "zcat $infile |" or die ("Could not open $infile for reading");
  }
  else {
    open(IN,'<',$infile) or die ("Could not open $infile for reading");
  }
  
  # Read through the file and parse out the desired fields
  while (<IN>) {
    chomp;
    
    next if /^MarkerName/;
    my @split = split(/\s+/);
    
    # Skip the risk allele if the variant is "0000"
    my $data = {
      'id'              => $split[0],
      'risk_allele'     => uc($split[1]),
      'p_value'         => $split[4],
      'description'     => $global_phenotype_description,
    };
    
    push(@phenotypes,$data);
  }
  
  close(IN);
  
  my %result = ('phenotypes' => \@phenotypes);
  return \%result;
}

sub parse_magic {
  my $infile = shift;

  my @phenotypes;

  # Open the input file for reading
  if($infile =~ /gz$/) {
    open IN, "zcat $infile |" or die ("Could not open $infile for reading");
  }
  else {
    open(IN,'<',$infile) or die ("Could not open $infile for reading");
  }

  my @headers;

  # Read through the file and parse out the desired fields
  while (<IN>) {
    chomp;

    if(/^snp/i) {
      @headers = split /\s+/;
    }
    else {
      $DB::single = 1;

      my @split = split(/\s+/);

      my %hash = map {$headers[$_] => $split[$_]} 0..$#split;

      # fix p-value
      $hash{pvalue} =~ s/^\./0\./;

      my $data = {
        'id'              => $hash{snp} || $hash{Snp},
        'risk_allele'     => uc($hash{effect_allele}),
        'p_value'         => $hash{pvalue},
        'description'     => $global_phenotype_description,
      };

      push(@phenotypes,$data);
    }
  }

  close(IN);

  my %result = ('phenotypes' => \@phenotypes);
  return \%result;
}

sub parse_orphanet {
  my $infile = shift;
  my $core_db_adaptor = shift;
  
  my $ga = $core_db_adaptor->get_GeneAdaptor;
  die("ERROR: Could not get gene adaptor") unless defined($ga);
  
  my @phenotypes;
  
  my %orphanet_special_char = (
    '\\\E7' => 'c',
    '\\\C5' => 'A',
    '\\\E8' => 'e',
    '\\\E9' => 'e',
    '\\\EB' => 'e',
    '\\\F6' => 'o',
    '\\\FC' => 'u',
    '\\\ED' => 'i',
    '\\\E4' => 'a',
  );

  my $xml_parser   = XML::LibXML->new();
  my $orphanet_doc = $xml_parser->parse_file($infile);
  
  foreach my $disorder ($orphanet_doc->findnodes('JDBOR/DisorderList/Disorder')) {
    my ($orpha_number_node) = $disorder->findnodes('./OrphaNumber');
    my $orpha_number = $orpha_number_node->to_literal;
    my ($name_node) = $disorder->findnodes('./Name');
    my $name = $name_node->to_literal;

    # Replace special characters
    if ($name =~ /\\/) {
      foreach my $char (keys(%orphanet_special_char)) {
        my $new_char = $special_characters{$char};
        $name =~ s/$char/$new_char/g;
      }
    }

    my @gene_nodes = $disorder->findnodes('./DisorderGeneAssociationList/DisorderGeneAssociation/Gene');
    
    foreach my $gene_node(@gene_nodes) {
      my $ref;
      my %ens_ids;
      
      #get the HGNC xref
      foreach my $external_reference_node ($gene_node->findnodes('./ExternalReferenceList/ExternalReference')) {
        my ($source_node) = $external_reference_node->findnodes('./Source');
        if ($source_node->to_literal =~ /HGNC/) {
          my ($ref_node) = $external_reference_node->findnodes('./Reference');
          $ref = $ref_node->to_literal;
        }
        if ($source_node->to_literal =~ /ensembl/i) {
          my ($ens_node) = $external_reference_node->findnodes('./Reference');
          $ens_ids{$ens_node->to_literal} = 1;
        }
      }
      if (defined($ref) || scalar keys %ens_ids) {
        my $genes = [];
        
        if(scalar keys %ens_ids) {
          foreach my $ens_id(keys %ens_ids) {
            my $g = $ga->fetch_by_stable_id($ens_id);
            push @$genes, $g if defined($g);
          }
        }

        if(scalar @$genes == 0 && defined($ref)) {
          my $genes = $ref =~ /^\d+$/ ? $ga->fetch_all_by_description('%HGNC:'.$ref.']%') : $ga->fetch_all_by_external_name($ref, 'HGNC');
          
          # we don't want any LRG genes
          @$genes = grep {$_->stable_id !~ /^LRG_/} @$genes;
          
          # and we don't want any duplicates if possible
          @$genes = grep {$ens_ids{$_->stable_id}} @$genes if scalar keys %ens_ids;
        }
        
        if ( scalar(@$genes) > 0) {
          if(scalar @$genes != 1) {
            print STDERR "WARNING: Found ".(scalar @$genes)." matching Ensembl genes for HGNC ID $ref or Ens IDs ".join(",", keys %ens_ids)."\n";
          }
          
          next unless scalar @$genes;
          
          foreach my $gene(@$genes) {
            push @phenotypes, {
              'id' => $gene->stable_id,
              'description' => $name,
              'external_id' => $orpha_number,
              'accessions'  => [ 'Orphanet:' . $orpha_number],
               ontology_mapping_type =>'is', 
              'seq_region_id' => $gene->slice->get_seq_region_id,
              'seq_region_start' => $gene->seq_region_start,
              'seq_region_end' => $gene->seq_region_end,
              'seq_region_strand' => $gene->seq_region_strand
            };
          }
        }
      }
    }
  }
  
  return {'phenotypes' => \@phenotypes};
}

sub parse_ddg2p {
  my $infile = shift;
  my $core_db_adaptor = shift;
  
  my $ga = $core_db_adaptor->get_GeneAdaptor;
  die("ERROR: Could not get gene adaptor") unless defined($ga);
  
  my @phenotypes;
  my $fh;
  
  # Open the input file for reading
  if($infile =~ /gz$/) {
    open $fh, "zcat $infile |" or die ("Could not open $infile for reading");
  }
  else {
    open($fh,'<',$infile) or die ("Could not open $infile for reading");
  }
  
  my %headers;
  
  my $csv = Text::CSV->new ( { binary => 1 } )  # should set binary attribute.
            or die "Cannot use CSV: ".Text::CSV->error_diag ();

  # Get columns headers
  $csv->column_names ($csv->getline ($fh));

  # Read through the file and parse out the desired fields
  while (my $content = $csv->getline_hr ($fh)) {
      
    # get data from the line
    my $symbol  = $content->{'gene symbol'};
    my $allelic = $content->{'allelic requirement'};
    my $mode    = $content->{'mutation consequence'};
    my $phen    = $content->{'disease name'};
    my $id      = $content->{'disease mim'};
    my @accns   = split/\;/,$content->{'phenotypes'};

    if ($symbol && $phen) {
      $phen =~ s/\_/ /g;

      my $genes = $ga->fetch_all_by_external_name($symbol, 'HGNC');

      # we don't want any LRG genes
      @$genes = grep {$_->stable_id !~ /^LRG_/} @$genes;

      # try restricting by name
      if (scalar(@$genes) > 1) {
        my @tmp = grep {$_->external_name eq $symbol} @$genes;
        $genes = \@tmp if scalar @tmp;
      }

      if (scalar @$genes != 1) {
        $DB::single = 1;
        print STDERR "WARNING: Found ".(scalar @$genes)." matching Ensembl genes for HGNC ID $symbol\n";
      }

      next unless scalar @$genes;

      foreach my $gene(@$genes) {
        push @phenotypes, {
          'id' => $gene->stable_id,
          'description' => $phen,
          'external_id' => $id,
          'seq_region_id' => $gene->slice->get_seq_region_id,
          'seq_region_start' => $gene->seq_region_start,
          'seq_region_end' => $gene->seq_region_end,
          'seq_region_strand' => $gene->seq_region_strand,
          'mutation_consequence' => $mode,
          'inheritance_type' => $allelic,
          'accessions' => \@accns,
          ontology_mapping_type =>'involves' 
        };
      }
    }
  }
  close($fh);

  return {'phenotypes' => \@phenotypes};
}

sub parse_cancer_gene_census {
  my $infile = shift;
  my $core_db_adaptor = shift;

  my $ga = $core_db_adaptor->get_GeneAdaptor;
  die("ERROR: Could not get gene adaptor") unless defined($ga);

  my @phenotypes;

  # Open the input file for reading
  if($infile =~ /gz$/) {
    open IN, "zcat $infile |" or die ("Could not open $infile for reading");
  }
  else {
    open(IN,'<',$infile) or die ("Could not open $infile for reading");
  }

  # Read through the file and parse out the desired fields
  while (<IN>) {
    chomp;

    my @row_data = split(/\t/, $_);

    # get data
    my $gene_id   = $row_data[1];
    my $phen      = $row_data[4];
    my $accession = $row_data[5];
    my $pmids     = $row_data[6];

    ## change accession format
    $accession =~ s/\_/\:/;

    my $gene = $ga->fetch_by_stable_id($gene_id);

    if(!$gene) {
      print STDERR "WARNING: Ensembl gene $gene_id not found in the Ensembl Core database. Skipped!\n";
      next;
    }

    my %data = (
      'id' => $gene_id,
      'description' => $phen,
      'seq_region_id' => $gene->slice->get_seq_region_id,
      'seq_region_start' => $gene->seq_region_start,
      'seq_region_end' => $gene->seq_region_end,
      'seq_region_strand' => $gene->seq_region_strand,
      'accessions' => [$accession],
      'ontology_mapping_type' => 'is'
    );

    $data{'study'} = $pmids if ($pmids);

    push(@phenotypes,\%data);
  }
  close(IN);

  return {'phenotypes' => \@phenotypes};
}


sub parse_omim_gene {
  my $infile = shift;
  open(IN, ($infile =~ /(z|gz)$/i ? "zcat $infile | " : $infile)) or die ("Could not open $infile for reading");
  
  # first record is empty
  <IN>;
  
  # Read through the file and parse out the desired fields
  while (<IN>) {
    chomp;
    
    my @content = split(/\t/,$_);
    
    my $desc = $content[6];
    next if (!$desc || $desc eq '');
    $desc =~ s/^\s+//;
    $desc = (split(';', $desc))[0]; # Only keep the first term

    my $gene = $content[0];
    next if (!$gene || $gene eq '');
    
    push @phenotypes, {
      'id'                => $gene,
      'description'       => $desc,
      'external_id'       => $content[5],
      'seq_region_id'     => $content[1],
      'seq_region_start'  => $content[2],
      'seq_region_end'    => $content[3],
      'seq_region_strand' => $content[4],
      'type'              => 'Gene',
    };
  }
  close IN;
  
  return {'phenotypes' => \@phenotypes}; 
}
 
#sub parse_omim_gene {
#  my $infile = shift;
#  my $mim2gene_file = shift;
#  my $core_db_adaptor = shift;
#  
#  # parse mim2gene
#  my $mim2gene = parse_mim2gene($mim2gene_file);
#  
#  my $ga = $core_db_adaptor->get_GeneAdaptor;
#  die("ERROR: Could not get gene adaptor") unless defined($ga);
# 
#  my @phenotypes;
#  
#  # set input record separator
#  local $/ = "*RECORD*";
#  
#  # Open the input file for reading
#  open(IN, ($infile =~ /(z|gz)$/i ? "zcat $infile | " : $infile)) or die ("Could not open $infile for reading");
#  
#  # first record is empty
#  <IN>;
#  
#  # Read through the file and parse out the desired fields
#  while (<IN>) {
#    chomp;
#    
#    if(/\*FIELD\*\s+NO\n(\d+)/){
#      my $number = $1;
#      next unless $mim2gene->{$number} && $mim2gene->{$number}->{symbols} ne '-';
#      
#      if(/\*FIELD\*\sTI\n([\^\#\%\+\*]*)\d+(.*)\n/){
#        my $label = $2; # taken from description as acc is meaning less
#        my $type = $1;
#        $label =~ s/\;\s[A-Z0-9]+$//; # strip gene name at end
#        $label = $label." [".$type.$number."]";
#        
#        my $symbol = $mim2gene->{$number}->{symbols};
#        my $genes = $ga->fetch_all_by_external_name($symbol, 'HGNC');
#        
#        # we don't want any LRG genes
#        @$genes = grep {$_->stable_id !~ /^LRG_/} @$genes;
#        
#        # try to get only the one with the correct display label
#        if(scalar @{$genes} > 0) {
#          my @tmp = grep {$_->external_name eq $symbol} @$genes;
#          @$genes = @tmp if scalar @tmp;
#        }
#        
#        if ( scalar(@$genes) > 0) {
#          if(scalar @$genes != 1) {
#            print STDERR "WARNING: Found ".(scalar @$genes)." matching Ensembl genes for HGNC ID $symbol\n";
#          }
#          
#          next unless scalar @$genes;
#          
#          foreach my $gene(@$genes) {
#            push @phenotypes, {
#              'id' => $gene->stable_id,
#              'description' => $label,
#              'external_id' => $number,
#              'seq_region_id' => $gene->slice->get_seq_region_id,
#              'seq_region_start' => $gene->seq_region_start,
#              'seq_region_end' => $gene->seq_region_end,
#              'seq_region_strand' => $gene->seq_region_strand,
#              'type' => 'Gene',
#            };
#          }
#        }
#      }
#    }
#    
#    last if scalar @phenotypes >= 20;
#  }
#  
#  return {'phenotypes' => \@phenotypes};
#}

sub parse_mim2gene {
  my $file = shift;
  
  # Open the input file for reading
  if($infile =~ /gz$/) {
    open IN, "zcat $infile |" or die ("Could not open $infile for reading");
  }
  else {
    open(IN,'<',$infile) or die ("Could not open $infile for reading");
  }
  
  my %return;
  
  while(<IN>) {
    next if /^\#/;
    chomp;
    
    my ($mim, $type, $genes, $symbols) = split /\s+/;
    
    $return{$mim} = {
      type => $type,
      genes => $genes,
      symbols => $symbols
    };
  }
  
  close IN;
  
  return \%return;
}

sub parse_mim_dump {
  my $file = shift;
  
  # Open the input file for reading
  if($infile =~ /gz$/) {
    open IN, "zcat $infile |" or die ("Could not open $infile for reading");
  }
  else {
    open(IN,'<',$infile) or die ("Could not open $infile for reading");
  }
  
  my @phenotypes;
  
  while(<IN>) {
    next if /^stable_id/;
    chomp;
    
    my ($gene, $sr, $s, $e, $str, $id, $desc) = split /\t/;
    
    next unless $desc && $gene;
    
    # sometimes they have spaces at the start?
    $desc =~ s/^\s+//;
    
    push @phenotypes, {
      'id' => $gene,
      'description' => $desc,
      'external_id' => $id,
      'seq_region_id' => $sr,
      'seq_region_start' => $s,
      'seq_region_end' => $e,
      'seq_region_strand' => $str,
      'type' => 'Gene',
    }
  }
  
  return {'phenotypes' => \@phenotypes};
}

# dbGaP data
sub parse_dbgap {
  my $infile = shift;
  
  my @phenotypes;
  my $study_prefix = 'pha';
  my @study_0 = ('000000','00000','0000','000','00','0',''); # Fill the study ID up to 6 numbers after the prefix "pha"
  
  # Open the input file for reading
  open(IN,'<',$infile) or die ("Could not open $infile for reading");
  
  # Read through the file and parse out the desired fields
  while (<IN>) {
    chomp;
    
    my @content = split(/\t/,$_);
    
    next if ($content[0] !~ /^dbgap$/i);
    
    my $phenotype = $content[1];
    my $rs_id     = $content[2];
    my $pvalue    = $content[3];
    my $gene      = $content[6];
    my $pubmed_id = ($content[8] ne '') ? $content[8] : undef;
    my $study     = ($content[21] ne '') ? $content[21] : undef;
       $study = $study_prefix.$study_0[length($study)].$study if (defined($study));

    $gene =~ s/ /,/g;

    # Store 2 decimals
    if ($pvalue =~ /^(\d\.\d+)(e-\d+)$/) {
      my $val = sprintf("%.2f", $1);
      $pvalue = "$val$2";
    }

    my %data = (
      'description'       => $phenotype,
      'associated_gene'   => $gene,
      'p_value'           => $pvalue,
      'study_description' => $study
    );
    
    # Parse the ids
    my @ids;
    $rs_id ||= "";
    while ($rs_id =~ m/(rs[0-9]+)/g) {
      push(@ids,$1);
    }
    $data{'variation_names'} = join(',',@ids);
    $data{'study'} = $pubmed_prefix . $pubmed_id if (defined($pubmed_id));
    
    # If we didn't get any rsIds, skip this row (this will also get rid of the header)
    warn("Could not parse any rsIds from string '$rs_id'") if (!scalar(@ids));
    next if (!scalar(@ids));
    
    map {
      my %t_data = %{\%data};
      $t_data{'id'} = $_;
      push(@phenotypes,\%t_data)
    } @ids;
  }
  close(IN);
  
  my %result = ('phenotypes' => \@phenotypes);
  
  return \%result;
}

sub parse_zfin {
  my $infile = shift;
  my $core_db_adaptor = shift;
  
  my $ga = $core_db_adaptor->get_GeneAdaptor;
  die("ERROR: Could not get gene adaptor") unless defined($ga);
  
  my @phenotypes;
  
  # Open the input file for reading
  if($infile =~ /gz$/) {
    open IN, "zcat $infile |" or die ("Could not open $infile for reading");
  }
  else {
    open(IN,'<',$infile) or die ("Could not open $infile for reading");
  }
  
  # Read through the file and parse out the desired fields
  while (<IN>) {
    chomp;
    
    my @data = split /\t/, $_;
    
    my $symbol  = $data[1];
    my $gene_id = $data[2];
    my $phen    = $data[4];
    
    for my $i(6, 8, 10, 12) {
      $phen .= ($phen ? ', ' : '').$data[$i] if $data[$i];
    }
    
    if($symbol && $phen) {
      my $genes = $ga->fetch_all_by_external_name($symbol);
      
      # try restricting by name
      if(scalar(@$genes) > 1) {
        my @tmp = grep {$_->external_name eq $symbol} @$genes;
        $genes = \@tmp if scalar @tmp;
      }
      
      if(scalar @$genes != 1) {
        print STDERR "WARNING: Found ".(scalar @$genes)." matching Ensembl genes for gene ID $symbol\n";
      }
      
      next unless scalar @$genes;
      
      foreach my $gene(@$genes) {
        push @phenotypes, {
          'id' => $gene->stable_id,
          'description' => $phen,
          'external_id' => $gene_id,
          'seq_region_id' => $gene->slice->get_seq_region_id,
          'seq_region_start' => $gene->seq_region_start,
          'seq_region_end' => $gene->seq_region_end,
          'seq_region_strand' => $gene->seq_region_strand,
        };
      }
    }
  }
  
  return {'phenotypes' => \@phenotypes};
}

sub get_mouse_phenotype_data_source_ids {
  my $phenotype_file = shift;
  my $data_source = shift;
  my $version = shift; 
  my $dbh = $db_adaptor->dbc->db_handle;
  my $fh = FileHandle->new($phenotype_file, 'r');
  my $source_names = {};
  my $source_name2id = {};

  while (<$fh>) {
    chomp;
    my @pairs = split("\t", $_);
    my $hash;
    foreach my $pair (@pairs) {
      my ($key, $value) = split('=', $pair);
      if ($key eq 'resource_name' && $data_source eq 'impc') {
        $source_names->{$value} = 1;
      }
      if ($key eq 'project_name' && $data_source eq 'mgi') {
        $source_names->{$value} = 1;
      }
    }
  }
  $fh->close();
  foreach my $source_name (keys %$source_names) {
    my $stmt = qq{ SELECT source_id FROM source WHERE name = '$source_name' LIMIT 1};
    my $sth = $dbh->prepare($stmt);
    $sth->execute();
    my $source_id;
    $sth->bind_columns(\$source_id);
    $sth->fetch();
    if (defined($source_id)) {
      $source_name2id->{$source_name} = $source_id;
    } else {
      warn "Could not fetch dbID for source name $source_name\n";
    }
    $sth->finish();
  }

  return $source_name2id;
}

sub clear_mouse_phenotype_data_from_last_release {
  my $source_name2id = shift;
  my $dbh = $db_adaptor->dbc->db_handle;
  my $source_ids = join(',', values %$source_name2id);
  $dbh->do(qq{ DELETE pfa FROM phenotype_feature_attrib pfa JOIN phenotype_feature pf ON pfa.phenotype_feature_id = pf.phenotype_feature_id AND pf.source_id IN ($source_ids);} );
  $dbh->do(qq{ DELETE FROM phenotype_feature WHERE source_id IN ($source_ids);} );
}

sub parse_mouse_phenotype_data {
  my $infile = shift;
  my $marker_coords = shift;
  my $data_source = shift;
  my $source_name2ids = shift;

  my $individual_adaptor = $db_adaptor->get_IndividualAdaptor;

  my @phenotypes = ();
  my $already_inserted = {};

  my @data_attrib = (
    'allele_symbol',
    'allele_accession_id',
    'marker_accession_id',);

  my $fh = FileHandle->new($infile, 'r');
  while (<$fh>) {
    chomp;
    my %data = ();
    my @pairs = split("\t", $_);
    my $hash;
    foreach my $pair (@pairs) {
      my ($key, $value) = split("=", $pair);
      if ($key eq 'allele_symbol') { # change e.g. Cdk5rap2<tm1a(EUCOMM)Wtsi> to Cdk5rap2_tm1a(EUCOMM)Wtsi
        $value =~ s/</_/g;
        $value =~ s/>//g;
      }
      $hash->{$key} = $value;
    }

    foreach my $value (@data_attrib) {
      $data{$value} = $hash->{$value};
    }
    my $marker_accession_id = $hash->{marker_accession_id};
    my $description = $hash->{mp_term_name};
    if (!$description) {
      warn "No description/mp_term_name for $_\n";
      next;
    }
    $data{description} = $description;
    $data{accession} = $hash->{mp_term_id};

    my $source = '';
    if ($data_source eq 'impc') {
      $source = $hash->{resource_name};
    } else {
      $source = $hash->{project_name};
    }
    $data{source} = $source;
    $data{source_id} = $source_name2ids->{$source};

    my $pf_data           = $marker_coords->{$marker_accession_id};
    my $type              = $pf_data->{type};
    my $seq_region_id     = $pf_data->{seq_region_id};
    my $seq_region_start  = $pf_data->{seq_region_start};
    my $seq_region_end    = $pf_data->{seq_region_end};
    my $seq_region_strand = $pf_data->{seq_region_strand};

    $data{seq_region_id} = $seq_region_id;
    $data{seq_region_start} = $seq_region_start;
    $data{seq_region_end} = $seq_region_end;
    $data{seq_region_strand} = $seq_region_strand; 

    my $strain_name       = $hash->{strain_name};
    my $gender            = $hash->{sex};
    if ($gender) {
      unless ($gender eq 'female' || $gender eq 'male') {
        $gender = 'unknown';
      }
    }

    my $strain_id = 0;
    my $create_new_individual = 0;
    if ($strain_name) {
      # attribs
      my $p_value;
      my $strains = $individual_adaptor->fetch_all_by_name($strain_name);
      if (scalar @$strains > 0) {
        my @unique = ();
        foreach my $strain (@$strains) {
          if (lc $strain->gender eq $gender) {
            push @unique, $strain;
          }
        }
        if (scalar @unique > 0) {
          $strain_id = $unique[0]->dbID();
        } else {
          $create_new_individual = 1;
        }
      }
      if ($create_new_individual) {
        my $individual = new Bio::EnsEMBL::Variation::Individual(
            -name => $strain_name,
            -gender => $gender,
            -individual_type_id => 1,
        );
        $individual_adaptor->store($individual);
        $strain_id = $individual_adaptor->last_insert_id();
      }
    }
    $data{associated_gene} = $hash->{marker_symbol};

    if ($hash->{p_value}) {
      $data{p_value} = convert_p_value($hash->{p_value});
    }
    $data{external_id} = $hash->{mp_term_id};

    if ($strain_id) {
      $data{strain_id} = $strain_id;
    }

    my @object_ids = ();
    if ($pf_data->{object_id}) {
      @object_ids = split(';', $pf_data->{object_id});
      foreach my $object_id (@object_ids) {
        my $accession = $data{accession};
        my $source = $data{source};
        my $strain_id = $data{strain_id};
        my $key = join('-', ($accession, $source, $object_id, $strain_id));
        unless ($already_inserted->{$key}) {
          my %phenotype_data = ();
          foreach my $key (keys %data) {
            $phenotype_data{$key} = $data{$key};
          }
          $phenotype_data{id} = $object_id;   
          $phenotype_data{accessions} = [$accession];
          $phenotype_data{ontology_mapping_type} = 'is'; 
          push(@phenotypes, \%phenotype_data);
        }
        $already_inserted->{$key} = 1;
      }
    }
  }
  my %result = ('phenotypes' => \@phenotypes);
  return \%result;
}

sub get_attrib_types {
  my $db_adaptor = shift;
  
  my $sth = $db_adaptor->dbc->prepare(qq{
    SELECT code
    FROM attrib_type
  });
  $sth->execute();
  
  my ($attrib_type, @tmp_types);
  $sth->bind_columns(\$attrib_type);
  while ($sth->fetch()) {
    push @tmp_types, $attrib_type if ($attrib_type ne 'description');
  }
  $sth->finish;
  
  return \@tmp_types;
}

sub get_seq_region_ids {
  my $db_adaptor = shift;
  
  my $sth = $db_adaptor->dbc->prepare(qq{
    SELECT seq_region_id, name
    FROM seq_region
  });
  $sth->execute;
  
  my (%seq_region_ids, $id, $name);
  $sth->bind_columns(\$id, \$name);
  $seq_region_ids{$name} = $id while $sth->fetch();
  $sth->finish;
  
  return \%seq_region_ids;
}

sub get_dbIDs {
  my $rs_ids = shift;
  my $db_adaptor = shift;
  
  my $id_stmt = qq{
    SELECT
      DISTINCT
      v.variation_id,
      v.name
    FROM
      variation v,
      variation_feature vf
      
    WHERE
      v.name = ? AND
      v.variation_id=vf.variation_id
    LIMIT 1
  };
  my $syn_stmt = qq{
    SELECT
      DISTINCT
      v.variation_id,
      v.name
    FROM
      variation_feature vf,
      variation_synonym vs JOIN
      variation v ON vs.variation_id = v.variation_id
    WHERE
      vs.name = ? AND
      v.variation_id=vf.variation_id
    LIMIT 1
  };
  my $id_sth = $db_adaptor->dbc->prepare($id_stmt);
  my $syn_sth = $db_adaptor->dbc->prepare($syn_stmt);
  
  my %mapping;
  
  foreach my $rs_id (@{$rs_ids}) {
    $id_sth->bind_param(1,$rs_id,SQL_VARCHAR);
    $id_sth->execute();
    my ($var_id,$var_name);
    $id_sth->bind_columns(\$var_id,\$var_name);
    $id_sth->fetch();
    
    # If we couldn't find the rs_id, look in the synonym table
    if (!defined($var_id)) {
      $syn_sth->bind_param(1,$rs_id,SQL_VARCHAR);
      $syn_sth->execute();
      $syn_sth->bind_columns(\$var_id,\$var_name);
      $syn_sth->fetch();
    }
    warn "$rs_id - no mapping found in variation db\n" unless $var_id ;  
    $mapping{$rs_id} = [$var_id,$var_name] if $var_id && $var_name;
  }
  
  return \%mapping;
}

sub get_coords {
  my $ids = shift;
  my $type = shift;
  my $db_adaptor = shift;
  
  my $tables;
  my $where_clause;
  
  my @object_ids = ($type eq 'Variation') ? map { $variation_ids->{$_}[1] } keys(%$variation_ids): @$ids;

  if($type eq 'Variation') {
    $tables = 'variation_feature f, variation v';
    $where_clause = 'v.variation_id = f.variation_id AND v.name = ?';
  }
  elsif($type =~ /Structural/) {
    $tables = 'structural_variation_feature f, structural_variation v';
    $where_clause = 'v.structural_variation_id = f.structural_variation_id AND v.variation_name = ?';
  }
  elsif($type eq 'Gene') {
    $tables = 'gene f';
    $where_clause = 'f.stable_id = ?';
  }
  
  my $sth = $db_adaptor->dbc->prepare(qq{
    SELECT
      f.seq_region_id,
      f.seq_region_start,
      f.seq_region_end,
      f.seq_region_strand
    FROM
      $tables
    WHERE
      $where_clause
  });
  
  my $coords = {};
  my ($sr_id, $start, $end, $strand);
  
  foreach my $id (@object_ids) {
    $sth->bind_param(1,$id,SQL_VARCHAR);
    $sth->execute();
    $sth->bind_columns(\$sr_id, \$start, \$end, \$strand);
    
    push @{$coords->{$id}}, {
      seq_region_id => $sr_id,
      seq_region_start => $start,
      seq_region_end => $end,
      seq_region_strand => $strand
    } while $sth->fetch();
  }
  
  $sth->finish();
  
  return $coords;
}

sub get_or_add_source {
  my $source_name = shift;
  my $source_description = shift;
  my $source_url    = shift;
  my $source_status = shift;
  my $db_adaptor    = shift;
  
  my $stmt = qq{
    SELECT
      source_id
    FROM
      source
    WHERE
      name = '$source_name'
    LIMIT 1
  };
  my $sth = $db_adaptor->dbc->prepare($stmt);
  $sth->execute();
  my $source_id;
  $sth->bind_columns(\$source_id);
  $sth->fetch();
  
  if (!defined($source_id)) {
    $stmt = qq{
      INSERT INTO
        source (
          name,
          description,
          url,
          somatic_status,
          version
        )
      VALUES (
        '$source_name',
        '$source_description',
        '$source_url',
        '$source_status',
        $source_version
      )
    };
    $db_adaptor->dbc->do($stmt);
    $source_id = $db_adaptor->dbc->db_handle->{'mysql_insertid'};
    
    print STDOUT "Added source for $source_name (source_id = $source_id)\n" if ($verbose);
  }
  else {
    $stmt = qq{
      UPDATE
        source 
      SET name=?,
        description=?,
        url=?,
        version=?
      WHERE
        source_id=?
    };
    my $update_source_sth =$db_adaptor->dbc->prepare($stmt);
    $update_source_sth->execute($source_name,$source_description,$source_url,$source_version,$source_id);
  }

  return $source_id;
}

sub add_phenotypes {
  my $phenotypes = shift;
  my $coords = shift;
  my $source_id = shift;
  my $object_type = shift;
  my $db_adaptor = shift;

  my $st_col = ($source =~ m/dbgap/i) ? 'name' : 'description';

  # Prepared statements
  my $st_ins_stmt = qq{
    INSERT INTO study (
      source_id,
      external_reference,
      study_type,
      $st_col
    )
    VALUES (
      $source_id,
      ?,
      ?,
      ?
    )
  };
  
  my $extra_cond = '';
  my $left_join = '';
  
  # add special joins for checking certain sources
  if ($source =~ m/uniprot/i) {
    $left_join = qq{
      LEFT JOIN
      (
        phenotype_feature_attrib pfa
        JOIN attrib_type at
        ON pfa.attrib_type_id = at.attrib_type_id
      )
      ON pf.phenotype_feature_id = pfa.phenotype_feature_id
    };
    
    $extra_cond = 'AND at.code = "variation_names" and pfa.value = ? ';
  }
  elsif ($source =~ m/omim/i) {
    $left_join = qq{
      LEFT JOIN
      (
        phenotype_feature_attrib pfa1
        JOIN attrib_type at1
        ON pfa1.attrib_type_id = at1.attrib_type_id
      )
      ON pf.phenotype_feature_id = pfa1.phenotype_feature_id
      
      LEFT JOIN
      (
        phenotype_feature_attrib pfa2
        JOIN attrib_type at2
        ON pfa2.attrib_type_id = at2.attrib_type_id
      )
      ON pf.phenotype_feature_id = pfa2.phenotype_feature_id
    };
    
    $extra_cond = qq{
      AND at1.code = "risk_allele"
      AND pfa1.value = ?
      AND at2.code = "associated_gene"
      AND pfa2.value = ?
    };
  }
  elsif ($source =~ m/nhgri/i) {
    $left_join = qq{
      LEFT JOIN
      (
        phenotype_feature_attrib pfa
        JOIN attrib_type at
        ON pfa.attrib_type_id = at.attrib_type_id
      )
      ON pf.phenotype_feature_id = pfa.phenotype_feature_id
    };
    
    $extra_cond = qq{ AND at.code = "p_value" AND pfa.value = ? };
  }
  
  my $pf_check_stmt = qq{
    SELECT
      pf.phenotype_feature_id
    FROM
      phenotype_feature pf
    $left_join
    WHERE
      pf.object_id = ? AND
      pf.type = ? AND
      pf.phenotype_id = ? AND
      pf.source_id = ? AND
      (pf.study_id = ? OR pf.study_id IS NULL)
      $extra_cond
    LIMIT 1
  };
  
  my $pf_ins_stmt = qq{
    INSERT INTO phenotype_feature (
      phenotype_id,
      source_id,
      study_id,
      type,
      object_id,
      is_significant,
      seq_region_id,
      seq_region_start,
      seq_region_end,
      seq_region_strand
    )
    VALUES (
      ?,
      ?,
      ?,
      ?,
      ?,
      ?,
      ?,
      ?,
      ?,
      ?
    )
  };
  
  my $attrib_ins_stmt = qq{
    INSERT INTO phenotype_feature_attrib (
      phenotype_feature_id,
      attrib_type_id,
      value
    )
    SELECT ?, at.attrib_type_id, ?
    FROM attrib_type at
    WHERE at.code = ?
  };
  
  my $attrib_ins_cast_stmt = qq{
    INSERT INTO phenotype_feature_attrib (
      phenotype_feature_id,
      attrib_type_id,
      value
    )
    SELECT ?, at.attrib_type_id, CAST(? AS CHAR)
    FROM attrib_type at
    WHERE at.code = ?
  };

  my $attrib_id_ext_stmt = qq{ 
    SELECT attrib_id from attrib, attrib_type
    WHERE attrib_type.code = 'ontology_mapping'
    AND attrib.attrib_type_id = attrib_type.attrib_type_id
    AND attrib.value = ?};
 
  my $ontology_accession_ins_stmt = qq{
    INSERT IGNORE INTO phenotype_ontology_accession
    (phenotype_id, accession, mapped_by_attrib, mapping_type)
    values (?,?,?,?)
   };
 
  my $st_ins_sth   = $db_adaptor->dbc->prepare($st_ins_stmt);
  my $pf_check_sth = $db_adaptor->dbc->prepare($pf_check_stmt);
  my $pf_ins_sth   = $db_adaptor->dbc->prepare($pf_ins_stmt);
  my $attrib_ins_sth = $db_adaptor->dbc->prepare($attrib_ins_stmt);
  my $attrib_ins_cast_sth = $db_adaptor->dbc->prepare($attrib_ins_cast_stmt);
  my $ontology_accession_ins_sth = $db_adaptor->dbc->prepare($ontology_accession_ins_stmt);

  ## get the attrib id for the type of description to ontology term linking
  my $attrib_id_ext_sth = $db_adaptor->dbc->prepare($attrib_id_ext_stmt);
  my $mapped_by = ($source eq 'Orphanet' ? 'Orphanet' : 'Data source');
  $attrib_id_ext_sth->execute($mapped_by);
  my $ont_attrib_type = $attrib_id_ext_sth->fetchall_arrayref();

  # First, sort the array according to the phenotype description
  my @sorted = sort {($a->{description} || $a->{name}) cmp ($b->{description} || $b->{name})} @{$phenotypes};
  my $study_count = 0;
  my $phenotype_feature_count = 0;
  
  my $total = scalar @sorted;
  my $i = 0;
  
  while (my $phenotype = shift(@sorted)) {
    progress($i++, $total);

    $object_type = $phenotype->{type} if defined($phenotype->{type});

    my $object_id = $phenotype->{"id"};
    # If the rs could not be mapped to a variation id, skip it
    next if $object_type =~ /Variation/ && (!defined($variation_ids->{$object_id}));

    $object_id = $variation_ids->{$object_id}[1] if ($object_type eq 'Variation');

    # if we have no coords, skip it
    my $study_id;
    
    if(defined($phenotype->{study}) || defined($phenotype->{study_description}) || defined($phenotype->{study_type})) {
      
      my $sql_study = '= ?';
      my $sql_type  = '= ?';
      
      # To avoid duplication of study entries
      if (!defined $phenotype->{"study"}) {$sql_study = 'IS NULL'; }
      if (!defined $phenotype->{"study_type"}) {$sql_type = 'IS NULL'; }  
      
      my $study_name = ($source =~ m/dbgap/i) ? ' AND name=? ' : '';
      
      my $st_check_stmt = qq{
        SELECT study_id
        FROM study
        WHERE
        source_id = $source_id AND
        external_reference $sql_study AND
        study_type $sql_type
        $study_name
        LIMIT 1
      };
         
      my $st_check_sth = $db_adaptor->dbc->prepare($st_check_stmt);
      my $param_num = 1;
      
      if (defined $phenotype->{"study"}) {
        $st_check_sth->bind_param($param_num,$phenotype->{"study"},SQL_VARCHAR);
        $param_num++;
      } 
      if (defined $phenotype->{"study_type"}) {
        $st_check_sth->bind_param($param_num,$phenotype->{"study_type"},SQL_VARCHAR);
        $param_num++;
      }
      $st_check_sth->bind_param($param_num,$phenotype->{"study_description"},SQL_VARCHAR) if ($source =~ m/dbgap/i);
      $st_check_sth->execute();
      $st_check_sth->bind_columns(\$study_id);
      $st_check_sth->fetch();
      
      if (!defined($study_id)) {
        $st_ins_sth->bind_param(1,$phenotype->{"study"},SQL_VARCHAR);
        $st_ins_sth->bind_param(2,$phenotype->{"study_type"},SQL_VARCHAR);
        $st_ins_sth->bind_param(3,$phenotype->{"study_description"},SQL_VARCHAR);
        $st_ins_sth->execute();
        $study_id = $db_adaptor->dbc->db_handle->{'mysql_insertid'};
        $study_count++;
      }
    }
    
    # Remove special characters from the phenotype description
    foreach my $char (keys(%special_characters)) {
      my $new_char = $special_characters{$char};
      $phenotype->{description} =~ s/$char/$new_char/g;
    }

    # get phenotype ID
    my $phenotype_id = get_phenotype_id($phenotype, $db_adaptor);

    foreach my $acc (@{$phenotype->{accessions}}){
      $acc =~ s/\s+//g;
      $acc = iri2acc($acc) if $acc =~ /^http/;
      $ontology_accession_ins_sth->execute( $phenotype_id, $acc,  $ont_attrib_type->[0]->[0], $phenotype->{ontology_mapping_type} ) ||die "Failed to import phenotype accession\n";
    }

    if ($phenotype->{"associated_gene"}) {
      $phenotype->{"associated_gene"} =~ s/\s//g;
      $phenotype->{"associated_gene"} =~ s/;/,/g;
    }

    $phenotype->{"p_value"} = convert_p_value($phenotype->{"p_value"}) if (defined($phenotype->{"p_value"}));

    # Check if this phenotype_feature already exists for this variation and source, in that case we probably want to skip it
    my $pf_id;

    $pf_check_sth->bind_param(1,$object_id,SQL_VARCHAR);
    $pf_check_sth->bind_param(2,$object_type,SQL_VARCHAR);
    $pf_check_sth->bind_param(3,$phenotype_id,SQL_INTEGER);
    $pf_check_sth->bind_param(4,$source_id,SQL_INTEGER);
    $pf_check_sth->bind_param(5,$study_id,SQL_INTEGER);
    # For uniprot data
    if ($source =~ m/uniprot/i) {
      $pf_check_sth->bind_param(6,$phenotype->{"variation_names"},SQL_VARCHAR);
    }
    # For nhgri-ebi gwas data
    elsif ($source =~ m/nhgri/i) {
      $pf_check_sth->bind_param(6,$phenotype->{"p_value"},SQL_VARCHAR);
    }
    # For omim data
    elsif ($source =~ m/omim/i) {
      $pf_check_sth->bind_param(6,$phenotype->{"risk_allele"},SQL_VARCHAR);
      $pf_check_sth->bind_param(7,$phenotype->{"associated_gene"},SQL_VARCHAR);
    }
    
    $pf_check_sth->execute();
    $pf_check_sth->bind_columns(\$pf_id);
    $pf_check_sth->fetch();
    next if (defined($pf_id));
    
    my $is_significant = defined($threshold) ? ($phenotype->{"p_value"} && $phenotype->{"p_value"} < $threshold ? 1 : 0) : 1;
    
    # Else, insert this variation annotation
    foreach my $coord(@{$coords->{$object_id}}) {
      $pf_ins_sth->bind_param(1,$phenotype_id,SQL_INTEGER);
      $pf_ins_sth->bind_param(2,$source_id,SQL_INTEGER);
      $pf_ins_sth->bind_param(3,$study_id,SQL_INTEGER);
      $pf_ins_sth->bind_param(4,$object_type,SQL_VARCHAR);
      $pf_ins_sth->bind_param(5,$object_id,SQL_VARCHAR);
      $pf_ins_sth->bind_param(6,$is_significant,SQL_INTEGER);
      $pf_ins_sth->bind_param(7,$coord->{seq_region_id},SQL_INTEGER);
      $pf_ins_sth->bind_param(8,$coord->{seq_region_start},SQL_INTEGER);
      $pf_ins_sth->bind_param(9,$coord->{seq_region_end},SQL_INTEGER);
      $pf_ins_sth->bind_param(10,$coord->{seq_region_strand},SQL_INTEGER);
      $pf_ins_sth->execute();
      $phenotype_feature_count++;
      
      # get inserted ID
      $pf_id = $db_adaptor->dbc->db_handle->{'mysql_insertid'};
      
      # add attribs
      foreach my $attrib_type(grep {defined($phenotype->{$_}) && $phenotype->{$_} ne ''} @attrib_types) {
        my $value = $phenotype->{$attrib_type};
        my $sth = $value =~ m/^\d+(\.\d+)?$/ ? $attrib_ins_cast_sth : $attrib_ins_sth;
        $sth->bind_param(1,$pf_id,SQL_INTEGER);
        $sth->bind_param(2,$value,SQL_VARCHAR);
        $sth->bind_param(3,$attrib_type,SQL_VARCHAR);
        $sth->execute();
      }
    }
  }
  end_progress();
  print STDOUT "$study_count new studies added\n" if ($verbose);
  print STDOUT "$phenotype_feature_count new phenotype_features added\n" if ($verbose);
}

sub get_phenotype_id {
  my $phenotype = shift;
  my $db_adaptor = shift;

  my ($name, $description);
  $name = $phenotype->{name};
  $description = $phenotype->{description};

  # Clean up
  $description =~ s/^\s+|\s+$//g; # Remove spaces at the beginning and the end of the description
  $description =~ s/\n//g; # Remove 'new line' characters
  
  die "ERROR: No description found for phenotype $name\n" if(!defined($description));

  # Check phenotype description in the format "description; name"
  if (!defined($name) || $name eq '') {
    my ($p_desc,$p_name) = split(";",$description);
    if ($p_name) {
      $p_name =~ s/ //g;
      if ($p_name =~ /^\w+$/) {
        $description = $p_desc;
        $name = $p_name;
      }
    }
  }
  
  if(scalar keys %phenotype_cache) {
    
    # check cache first
    return $phenotype_cache{$description} if defined $phenotype_cache{$description};
    
    my @tmp = keys %phenotype_cache;
    
    # lc everything
    my $description_bak = $description;
    $description = lc($description);
    
    # store mapped
    my %mapped;
    $mapped{lc($_)} = $_ for @tmp;
    @tmp = keys %mapped;
    
    # check if it matches lc
    if(defined($mapped{$description})) {
      $phenotype_cache{$description_bak} = $phenotype_cache{$mapped{$description}};
      return $phenotype_cache{$description_bak};
    }
    
    # try a fuzzy match using String::Approx
    my @matches = amatch($description, [10], @tmp);
  
    if(@matches) {
      
      # we only want the best match
      my $best = scalar @matches == 1 ? $matches[0] : (sort {abs(adist($description, $a)) <=> abs(adist($description, $b))} @matches)[0];
      
      # if distance is 0, return
      return $phenotype_cache{$mapped{$best}} if adist($description, $best) == 0;
      
      # find characters that differ
      my $diff = diff([split(//,$description)], [split(//,$best)]);
      
      my $skip = 0;
      
      # skip if mismatch is anything word-like
      foreach(map {@$_} @$diff) {
        $skip = 1 if $_->[2] =~ /\w/i;
      }
      
      # cache this match so we don't have to fuzz again
      $phenotype_cache{$description_bak} = $phenotype_cache{$mapped{$best}};
      
      return $phenotype_cache{$mapped{$best}} unless $skip;
    }
    
    # restore from backup before inserting
    $description = $description_bak;
  }


  # finally if no match, do an insert
  my $sth = $db_adaptor->dbc->prepare(qq{
    INSERT IGNORE INTO phenotype ( name, description ) VALUES ( ?,? )
  });
  
  $sth->bind_param(1,$name,SQL_VARCHAR);
  $sth->bind_param(2,$description,SQL_VARCHAR);
  $sth->execute();
  my $phenotype_id = $db_adaptor->dbc->db_handle->{'mysql_insertid'};
  
  # update cache
  $phenotype_cache{$description} = $phenotype_id;
  
  return $phenotype_id;
}

sub add_synonyms {
  my $synonyms = shift;
  my $variation_ids = shift;
  my $source_id = shift;
  my $db_adaptor = shift;
  
  # If we actually didn't get any synonyms, just return
  return if (!defined($synonyms) || !scalar(keys(%{$synonyms})));
  
  # Some prepeared statements needed for inserting the synonyms into database
  my $ins_stmt = qq{
    INSERT IGNORE INTO
      variation_synonym (
      variation_id,
      source_id,
      name
      )
    VALUES (
      ?,
      $source_id,
      ?
    )
  };
  my $ins_sth = $db_adaptor->dbc->prepare($ins_stmt);
  
  my $alt_count = 0;
  my $variation_count = 0;
  
  foreach my $rs_id (keys %{$variation_ids}) {
    
    my $var_id = $variation_ids->{$rs_id}[0];
    
    # If we have a variation id, we can proceed
    if (defined($var_id)) {
      
      $variation_count++;
      
      $ins_sth->bind_param(1,$var_id,SQL_INTEGER);
      
      # Handle all synonym ids for this rs_id
      while (my $alt_id = shift(@{$synonyms->{$rs_id}})) {
      
        # Add the id as synonym, if it is already present, it will just be ignored
        $ins_sth->bind_param(2,$alt_id,SQL_VARCHAR);
        $ins_sth->execute();
        $alt_count++;
      }
    }
  }
  
  print STDOUT "Added $alt_count synonyms for $variation_count rs-ids\n" if ($verbose);
}

sub add_set {
  my $set = shift;
  my $source_id = shift;
  my $db_adaptor = shift;
  my $set_from = shift;
  
  return if (!defined($set));
  return if (!defined($set_from) && ($set_from ne 'phenotype' || $set_from ne 'synonym')); 
   
  my $variation_set_id;
  
  # Get variation_set_id
  my $select_set_stmt = qq{
  SELECT v.variation_set_id
  FROM variation_set v, attrib a
  WHERE v.short_name_attrib_id=a.attrib_id 
  AND a.value = ?
  };
  my $sth1 = $db_adaptor->dbc->prepare($select_set_stmt);
  $sth1->bind_param(1,$set,SQL_VARCHAR);
  $sth1->execute();
  $sth1->bind_columns(\$variation_set_id);
  $sth1->fetch();
  return if (!defined($variation_set_id));
  
  # Insert into variation_set_variation
  my $insert_set_stmt = qq{
    INSERT IGNORE INTO variation_set_variation (variation_id,variation_set_id)
    SELECT distinct v.variation_id, ? 
  };
  if ($set_from eq 'phenotype') {
    $insert_set_stmt .= qq{ 
      FROM phenotype_feature pf, variation v WHERE 
        v.name=pf.object_id AND
        pf.type='Variation' AND
        pf.source_id=?
    };
  }
  elsif ($set_from eq 'synonym') {
    $insert_set_stmt .= qq{ 
      FROM variation_synonym vs, variation v WHERE 
        v.variation_id=vs.variation_id AND
        vs.source_id=?
    };
  }
  
  my $sth2 = $db_adaptor->dbc->prepare($insert_set_stmt);
  $sth2->bind_param(1,$variation_set_id,SQL_INTEGER);
  $sth2->bind_param(2,$source_id,SQL_INTEGER);
  $sth2->execute();
}


sub convert_p_value {

  my $pval = shift;
  
  my $sci_pval = '';
  # If a scientific format is not found, then ...
  if ($pval !~ /^\d+.*e.+$/i) {  
    # If a range format is found (e.g. 10^-2 > p > 10^-3)
    if ($pval =~ /^\d+\^(-\d+)/) {
      if (length("$1")==1) { $1 = "0$1"; } 
      $sci_pval = "1.00e$1"; # e.g 10^-2 > p > 10^-3 => 1.00e-2
    }
    # If a decimal format is found (e.g. 0.0023)
    elsif ($pval =~ /^\d+/){
      $sci_pval = $pval;
    #$sci_pval = sprintf("%.2e",$pval); # e.g. 0.002 => 2,30e-3
    }
    elsif ($pval =~ /^\w+/) {
      $sci_pval = "NULL";
    }
    
  }
  else {
    $pval =~ tr/E/e/;
    if ($pval =~ /^(\d+)(e-?\d+)/) {
      $pval="$1.00$2";  
    }
    if ($pval =~ /^(\d+\.\d{1})(e-?\d+)/) {
      $pval="$1"."0$2";  
    }
    $sci_pval = $pval;
  }
  return $sci_pval;
}


# Method to remove the carriage return character in each line of the input file
sub remove_carriage_return {
  my $infile = shift;
  my $replacement = shift;
  
  $replacement ||= "";
  
  my @ARGV_bak = @ARGV;
  @ARGV = ($infile);
  $^I = ".bak2";
  while (<>) {
    s/\r/$replacement/g;
    print;
  }
  # Restore the @ARGV variable
  @ARGV = @ARGV_bak;
}

# update or initiate progress bar
sub progress {
  my ($i, $total) = @_;
  
  my $width = 60;
  my $percent = int(($i/$total) * 100);
  my $numblobs = int((($i/$total) * $width) - 2);
  
  # this ensures we're not writing to the terminal too much
  return if defined($prev_prog) && $numblobs.'-'.$percent eq $prev_prog;
  $prev_prog = $numblobs.'-'.$percent;
  
  printf("\r% -${width}s% 1s% 10s", '['.('=' x $numblobs).($numblobs == $width - 2 ? '=' : '>'), ']', "[ " . $percent . "% ]");
}

## format ontology accessions
sub iri2acc{

  my $iri = shift;
  my @a = split/\//, $iri;
  my $acc = pop @a;
  $acc =~ s/\_/\:/;

  return $acc;
}

sub get_mouse_phenotype_data {
  my $working_dir = shift;
  my $data_source = shift;
  my $url = shift;
  my $phenotype_file = "$working_dir/$data_source\_phenotypes.txt";
  my $fh = FileHandle->new($phenotype_file, 'w');
  my $http = HTTP::Tiny->new();

  my $i      = 0;
  my $rows   = 0;
  my $server = 'http://www.ebi.ac.uk';
  my $ext     = "$url/select?q=*:*&rows=$i&wt=json";

  # get total number of rows
  while ($i == $rows) {
    $i += 1000;
    $ext = "$url/select?q=*:*&rows=$i&wt=json";
    my $response = $http->get($server.$ext, {
      headers => { 'Content-type' => 'application/json' }
    });
    my $hash = decode_json($response->{content});
    $rows = scalar @{$hash->{response}->{docs}};
  }
  my $start = 0;
  while ($start <= $rows) {
    $ext = "$url/select?q=*:*&rows=100&wt=json&start=$start";
    my $response = $http->get($server.$ext, {
      headers => { 'Content-type' => 'application/json' }
    });

    die "Failed!\n" unless $response->{success};
    if (length $response->{content}) {
      my $hash = decode_json($response->{content});
      foreach my $result (@{$hash->{response}->{docs}}) {
        my @pairs = ();
        foreach my $key (keys %$result) {
          my $value = $result->{$key};
          $value ||= '\N';
          push @pairs, $key . "=" . $value;
        }
        print $fh join("\t", @pairs), "\n";
      }
    }
    $start += 100;
  }
  $fh->close();
  return $phenotype_file;
}

sub update_mouse_phenotype_data_version {
  my $source_name2ids = shift;
  my $http = HTTP::Tiny->new();

  my $url = 'http://www.mousephenotype.org/data/release.json';
  my $response = $http->get($url, {
    headers => { 'Content-type' => 'application/json' }
  });
  my $hash = decode_json($response->{content});
  my $data_release_date = $hash->{data_release_date};
  die "data_release_date not defined\n" unless $data_release_date;

 # e.g. 20 October 2015
  my ($date, $month, $year) = split(' ', $hash->{data_release_date}) ;

  my $months = {
    'January' => '01',
    'February' => '02',
    'March' => '03',
    'April' => '04',
    'May' => '05',
    'June' => '06',
    'July' => '07',
    'August' => '08',
    'September' => '09',
    'October' => '10',
    'November' => '11',
    'December' => '12',
  };

  my $month_number = $months->{$month};
  die "month_number not defined\n" unless $month_number;

  my $version = sprintf '%d%02d%02d', $year, $month_number, $date;

  my $dbh = $db_adaptor->dbc->db_handle;

  foreach my $source_name (keys %$source_name2ids) {
    $dbh->do(qq{UPDATE source SET version=$version WHERE name='$source_name';}) or die $dbh->errstr;
  }
  return $version;
}

sub get_marker_coords {
  my $phenotype_file = shift;
  my $coord_file = shift;

  my $gene_adaptor = $core_db_adaptor->get_GeneAdaptor;

  my $markers;
  my $column_headers;
  # filter for marker accession ids (MGI)
  my $fh = FileHandle->new($phenotype_file, 'r');
  while (<$fh>) {
    chomp;
    my @pairs = split("\t", $_);
    my $hash;
    foreach my $pair (@pairs) {
      my ($key, $value) = split("=", $pair);
      $hash->{$key} = $value;
      $column_headers->{$key} = 1;
    }
    $markers->{$hash->{marker_accession_id}} = 1;
  }
  $fh->close();

  $fh = FileHandle->new($coord_file, 'r');
  my $header;
  my $marker_coords = {};
  while (<$fh>) {
    chomp;
    # test column names are still the same
    if (/^1\./) {
      my @values = split("\t", $_);
      foreach my $value (@values) {
        my ($column_number, $column_name) = split(/\s/, $value, 2);
        $header->{$column_number} = $column_name;
      }
      unless ($header->{'1.'} eq 'MGI Marker Accession ID') { die 'Header in MGI_MRK_Coord.rpt has changed.'};
      unless ($header->{'6.'} eq 'Chromosome') { die 'Header in MGI_MRK_Coord.rpt has changed.'};
      unless ($header->{'7.'} eq 'Start Coordinate') { die 'Header in MGI_MRK_Coord.rpt has changed.'};
      unless ($header->{'8.'} eq 'End Coordinate') { die 'Header in MGI_MRK_Coord.rpt has changed.'};
      unless ($header->{'9.'} eq 'Strand') { die 'Header in MGI_MRK_Coord.rpt has changed.'};
    }
   my @values = split("\t", $_);
    my $marker_acc = $values[0];
    if ($markers->{$marker_acc}) {
      unless ($marker_coords->{$marker_acc}) {
        my $genes = $gene_adaptor->fetch_all_by_external_name($marker_acc);
        if (scalar @$genes != 1) {
          my $number_of_genes = scalar @$genes;
          print STDERR "WARNING: Found $number_of_genes matching Ensembl genes for gene $marker_acc\n";
          my $unique_seq_region_id = {};
          foreach my $gene (@$genes) {
            $unique_seq_region_id->{$gene->slice->get_seq_region_id} = 1;
          }
          print STDERR "WARNING: genes $number_of_genes are on different slices\n" unless (scalar keys %$unique_seq_region_id == 1);
        }
        next unless scalar @$genes;

        my $marker_type = $values[1];
        print STDERR "WARNING: type is not Gene (type = $marker_type)\n" unless ($marker_type eq 'Gene');

        my $chromosome = $values[5];
        my $start_coord = $values[6];
        my $end_coord = $values[7];
        my $strand = $values[8];

        $marker_coords->{$marker_acc}->{type} = $marker_type;
        $marker_coords->{$marker_acc}->{object_id} = join(";", map {$_->stable_id} @$genes);
        $marker_coords->{$marker_acc}->{chromosome} = $chromosome;
        $marker_coords->{$marker_acc}->{seq_region_id} = $genes->[0]->slice->get_seq_region_id;
        $marker_coords->{$marker_acc}->{seq_region_start} = $start_coord;
        $marker_coords->{$marker_acc}->{seq_region_end} = $end_coord;
        $marker_coords->{$marker_acc}->{seq_region_strand} = ($strand eq '+') ? 1 : -1;
      }
    }
  }

  $fh->close();
  return $marker_coords;
}

# end progress bar
sub end_progress {
  progress(1,1);
  print "\n";
}


sub usage {
  
  print qq{
  Usage: perl import_phenotype_data.pl [OPTION]
  
  Import phenotype data into a Variation database
  
  Options:
    
      -verbose           Progress information is printed
      -help              Print this message
      
      -skip_phenotypes   Skip the study, phenotype_feature, phenotype_feature_attrib and phenotype tables insertions.
      -skip_synonyms     Skip the variation_synonym table insertion.
      -skip_sets         Skip the variation_set_variation table insertion.
      
    Database credentials are specified on the command line
    
      -host      Variation database host name (Required)
      -dbname    Variation database name (Required)
      -user      Variation database user (Required)
      -pass      Variation database password (Required)
      -port      Variation database port (Default: 3306)
    
      -chost     Core database host name (Required for gene-type sources)
      -cdbname   Core database name (Required for gene-type sources)
      -cuser     Core database user (Required for gene-type sources)
      -cpass     Core database password (Required for gene-type sources)
      -cport     Core database port (Default: 3306)
      
    An input file must be specified. This file contains the data that will be imported, typically tab-delimited
    and obtained from the UniProt or NHGRI-EBI GWAS catalog. If a new source is required, a method for parsing the
    source format into a "standard" data structure can be added.
    
      -infile          Typically a tab-delimited file (Required)
      
    The source of the data must be specified so that the correct parser and source table entry will be used. Currently
    supported sources are 'uniprot' and 'nhgri'.
    
      -source             String indicating the source of the data (Required)
      -version           Numerical version of the source (Required)
  } . "\n";
  exit(0);
}
