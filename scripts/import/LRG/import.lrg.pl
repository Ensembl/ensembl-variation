=pod

SYNOPSIS

  Script to perform various actions, related to LRGs, on a Core database

DESCRIPTION

  This script can be used to:
    - Import a LRG into the Core database
    - Add xrefs to genes on a LRG and to Ensembl genes, linking them to LRG genes
    - Add gene_attribs to Ensembl genes, indicating that they are completely or partially overlapped by a LRG
    - Remove a LRG from the Core database
  
EXAMPLE
  
  Display help message:
    perl import.lrg.pl -help
    
  Import a LRG and add xrefs, will download XML record from website:
    perl import.lrg.pl -host ens-genomics1 -port 3306 -user ******** -pass ********** -dbname homo_sapiens_core_58_37c -lrg_id LRG_1 -import -xrefs
    
  Add gene_attribs for Ensembl genes overlapping a LRG:
    perl import.lrg.pl -host ens-genomics1 -port 3306 -user ******** -pass ********** -dbname homo_sapiens_core_58_37c -lrg_id LRG_1 -overlap
    
  Clean a LRG from the Core database:
    perl import.lrg.pl -host ens-genomics1 -port 3306 -user ******** -pass ********** -dbname homo_sapiens_core_58_37c -lrg_id LRG_1 -clean
    
=cut

#!perl -w

use strict;

use Getopt::Long;
use List::Util qw (min max);
use LRG::LRG;
use LRG::LRGImport;
use LRG::LRGMapping;

use Bio::EnsEMBL::Registry;

# Some constants
my $LRG_COORD_SYSTEM_NAME = q{lrg};
my $LRG_BIOTYPE = q{LRG_gene};
my $LRG_ANALYSIS_LOGIC_NAME = q{LRG_import};
my $LRG_ANALYSIS_DESCRIPTION = q{Data from LRG database};
my $LRG_ANALYSIS_DISPLAY_LABEL = q{LRG Genes};
my $LRG_ANALYSIS_WEB_DATA = qq{{'colour_key' => 'rna_[status]','caption' => 'LRG gene','label_key' => '[text_label] [display_label]','name' => 'LRG Genes','default' => {'MultiTop' => 'gene_label','contigviewbottom' => 'transcript_label','MultiBottom' => 'collapsed_label','contigviewtop' => 'gene_label','alignsliceviewbottom' => 'as_collapsed_label','cytoview' => 'gene_label'},'multi_caption' => 'LRG genes'}};
my $HGNC_EXTERNAL_DB_NAME = q{HGNC};
my $LRG_EXTERNAL_DB_NAME = q{LRG};
my $LRG_EXTERNAL_STATUS = q{KNOWN};
my $LRG_EXTERNAL_PRIORITY = 10;
my $LRG_EXTERNAL_DB_DISPLAY_NAME = q{Locus Reference Genomic};
my $LRG_EXTERNAL_DB_RELEASE = 1;
my $LRG_EXTERNAL_DB_ACC_LINKABLE = 1;
my $LRG_EXTERNAL_DB_LABEL_LINKABLE = 0;
my $LRG_EXTERNAL_TYPE = q{MISC};
my $LRG_ENSEMBL_DB_NAME = q{ENS_LRG};
my $LRG_ENSEMBL_STATUS = q{KNOWN};
my $LRG_ENSEMBL_PRIORITY = 10;
my $LRG_ENSEMBL_DB_DISPLAY_NAME = q{LRG display in Ensembl};
my $LRG_ENSEMBL_DB_RELEASE = 1;
my $LRG_ENSEMBL_DB_ACC_LINKABLE = 1;
my $LRG_ENSEMBL_DB_LABEL_LINKABLE = 0;
my $LRG_ENSEMBL_TYPE = q{MISC};
my $LRG_EXTERNAL_XML = q{ftp://ftp.ebi.ac.uk/pub/databases/lrgex/};

my $host;
my $port;
my $user;
my $pass;
my $help;
my $clean;
my $verbose;
my $overlap;
my @lrg_ids;
my $input_file;
my $import;
my $add_xrefs;
my $max_values;
my $revert;
my $verify;
my $coredb;
my $otherfeaturesdb;
my $cdnadb;
my $vegadb;
my $purge;

usage() if (!scalar(@ARGV));

# get options from command line
GetOptions(
  'host=s'		=> \$host,
  'port=i'		=> \$port,
  'core=s'		=> \$coredb,
  'user=s'		=> \$user,
  'pass=s'		=> \$pass,
  'help!' 		=> \$help,
  'verbose!' 		=> \$verbose,
  'clean!' 		=> \$clean,
  'purge!'		=> \$purge,
  'overlap!' 		=> \$overlap,
  'lrg_id=s' 		=> \@lrg_ids,
  'input_file=s' 	=> \$input_file,
  'import!' 		=> \$import,
  'xrefs!' 		=> \$add_xrefs,
  'max!' 		=> \$max_values,
  'revert!' 		=> \$revert,
  'verify!'		=> \$verify,
  'otherfeatures=s'	=> \$otherfeaturesdb,
  'cdna=s'		=> \$cdnadb,
  'vega=s'		=> \$vegadb,
);

usage() if (defined($help));

die("Database credentials (-host, -port, -core, -user) need to be specified!") unless (defined($host) && defined($port) && defined($coredb) && defined($user));

# If an input XML file was specified, this will override any specified lrg_id. So get the identifier from within the file
if ((defined($import) || defined($verify) || defined($add_xrefs)) && defined($input_file)) {

  die("ERROR: Input file $input_file does not exist!") unless(-e $input_file);
  
  # create an LRG object from input file
  print STDOUT localtime() . "\tCreating LRG object from input XML file $input_file\n" if ($verbose);
  my $lrg = LRG::LRG::newFromFile($input_file) or die("ERROR: Could not create LRG object from XML file!");
  
  # find the LRG ID
  my $lrg_name = $lrg->findNode("fixed_annotation/id")->content();
  print STDOUT localtime() . "\tLRG ID is $lrg_name\n" if ($verbose);
  
  # Set the lrg_id array
  @lrg_ids = ($lrg_name);
}

# Check that the LRG id is on the correct format
die("Supplied LRG id is not in the correct format ('LRG_NNN')") if (grep($_ !~ m/^LRG\_[0-9]+$/,@lrg_ids));

# If doing something requiring the XML file but without specified input XML file or LRG ids, get a listing of published LRGs available at the ftp site
if ((defined($import) || defined($clean) || defined($overlap) || defined($verify) || defined($add_xrefs)) && !scalar(@lrg_ids)) {
  
  print STDOUT localtime() . "\tNo input XML file and no LRG id specified, fetching a LRG listing from the LRG server\n" if ($verbose);
  my $result = LRGImport::fetch_remote_lrg_ids([$LRG_EXTERNAL_XML]);
  
  if ($result->{$LRG_EXTERNAL_XML}{'success'}) {
    my @available = @{$result->{$LRG_EXTERNAL_XML}{'lrg_id'}};
    print "The following LRGs are available for import from the LRG public ftp server:\n";
    print "\t" . join("\n\t",@available) . "\n";
    my %entered;
    print "Enter the LRG ids you want to import (enter for all), enter a blank line when finished\n";
    my $id = 1;
    while ($id) {
      print "\tLRG id: ";
      $id = <>;
      chomp($id);
      if (grep($_ eq $id,@available)) {
	$entered{$id} = 1;
      }
      elsif (length($id) > 0) {
	print "\tLRG identifier not recognized!\n";
      }
    }
    
    if (scalar(keys(%entered)) == 0) {
      @lrg_ids = @available;
    }
    else {
      @lrg_ids = keys(%entered);
    }
    
    print "Will process " . join(", ",@lrg_ids) . "\n";
  }
  else {
    die("Could not get LRG listing from external db. Server said: " . $result->{$LRG_EXTERNAL_XML}{'message'});
  }
}

# Connect to core database
print STDOUT localtime() . "\tGetting human core db adaptor\n" if ($verbose);
my $dbCore = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -host => $host,
  -user => $user,
  -pass => $pass,
  -port => $port,
  -dbname => $coredb
) or die("Could not get a database adaptor to $coredb on $host:$port");
print STDOUT localtime() . "\tConnected to $coredb on $host:$port\n" if ($verbose);

$LRGImport::dbCore = $dbCore;

# If specified, connect to otherfeatures and cdna database
my $dbOther;
my $dbcDNA;
my $dbVega;
if (defined($otherfeaturesdb)) {
  print STDOUT localtime() . "\tGetting db adaptor for $otherfeaturesdb\n" if ($verbose);
  $dbOther = new Bio::EnsEMBL::DBSQL::DBAdaptor(
    -host => $host,
    -user => $user,
    -pass => $pass,
    -port => $port,
    -dbname => $otherfeaturesdb
  ) or die("Could not get a database adaptor to $otherfeaturesdb on $host:$port");
  print STDOUT localtime() . "\tConnected to $otherfeaturesdb on $host:$port\n" if ($verbose);
}
if (defined($cdnadb)) {
  print STDOUT localtime() . "\tGetting db adaptor for $cdnadb\n" if ($verbose);
  $dbcDNA = new Bio::EnsEMBL::DBSQL::DBAdaptor(
    -host => $host,
    -user => $user,
    -pass => $pass,
    -port => $port,
    -dbname => $cdnadb
  ) or die("Could not get a database adaptor to $cdnadb on $host:$port");
  print STDOUT localtime() . "\tConnected to $cdnadb on $host:$port\n" if ($verbose);
}
if (defined($vegadb)) {
  print STDOUT localtime() . "\tGetting db adaptor for $vegadb\n" if ($verbose);
  $dbVega = new Bio::EnsEMBL::DBSQL::DBAdaptor(
    -host => $host,
    -user => $user,
    -pass => $pass,
    -port => $port,
    -dbname => $vegadb
  ) or die("Could not get a database adaptor to $vegadb on $host:$port");
  print STDOUT localtime() . "\tConnected to $vegadb on $host:$port\n" if ($verbose);
}

# Put the db adaptors into an array
my @db_adaptors = ($dbCore,$dbOther,$dbcDNA,$dbVega);

#ÊGet a slice adaptor
print STDOUT localtime() . "\tGetting slice adaptor\n" if ($verbose);
my $sa = $dbCore->get_SliceAdaptor();

# Get the maximum key field values if required and print them
if ($max_values) {
  my $max_values = LRGImport::get_max_key();
  while (my ($field,$value) = each(%{$max_values})) {
    my ($table,$fld) = split(/\./,$field);
    print $table . "\t" . $fld . "\t" . $value . "\n";
  }
}

die("A tab-separated input file with table, field and max_value columns must be specified in order to revert the core db!") if (defined($revert) && !defined($input_file));
# Revert the database tables by deleting all rows with the specified field value above the specified maximum
if ($revert) {
  open(MV,'<',$input_file);
  while (<MV>) {
    chomp;
    my ($table,$field,$max_value) = split();
    LRGImport::remove_row([qq{$field > $max_value}],[$table]);
  }
  close(MV);
}

#ÊIf doing an import, check that the tables affected for adding the mapping information are sync'd across the relevant databases
if ($import) {
  my %max_increment = (
  'coord_system' 	=> 0,
  'meta'		=> 0,
  'seq_region'		=> 0
  );
  foreach my $table (keys(%max_increment)) {
  
    print STDOUT localtime() . "\tChecking Auto_increment value for $table\n" if ($verbose);
    #ÊCheck each db adaptor
    my $do_sync = 0;
    foreach my $dba (@db_adaptors) {
      # Skip this db adaptor if it's not defined
      next unless(defined($dba));
      
      my $stmt = qq{
        SHOW TABLE STATUS LIKE '$table'
      };
      my @current = keys(%{$dba->dbc->db_handle->selectall_hashref($stmt,'Auto_increment')}) or die("Could not get AUTO_INCREMENT value for " . $dba->dbc()->dbname() . ".$table");
      if ($max_increment{$table} > 0 && $max_increment{$table} != $current[0]) {
        $do_sync = 1;
      }
      print STDOUT localtime() . "\t\tAuto_increment value for " . $dba->dbc()->dbname() . ".$table is " . $current[0] . ($do_sync ? " => Needs to be sync'd!" : "") . "\n" if ($verbose);
      $max_increment{$table} = max($max_increment{$table},$current[0]);
    }
    # If the table needs to be sync'd, loop over the db adaptors and do this
    next unless($do_sync);
    
    foreach my $dba (@db_adaptors) {
      # Skip this db adaptor if it's not defined
      next unless(defined($dba));
      
      print STDOUT localtime() . "\tSyncing " . $dba->dbc()->dbname() . ".$table\n" if ($verbose);
      
      my $val = $max_increment{$table};
      my $stmt = qq{
        ALTER TABLE
          $table
        AUTO_INCREMENT = $val
      };
      $dba->dbc()->do($stmt) or die("Could not set AUTO_INCREMENT value for " . $dba->dbc()->dbname() . ".$table");
    }
  }
}
  
#ÊLoop over the specified LRG identifiers and process each one
while (my $lrg_id = shift(@lrg_ids)) {
  
  print localtime() . "\tProcessing $lrg_id\n";
  
  # Clean up data in the database if required
  if ($clean) {
    # Loop over all db adaptors
    foreach my $dba (@db_adaptors) {
      # Skip this db adaptor if it's not defined
      next unless(defined($dba));
	
      # Set the dbCore variable in LRGImport to the current db adaptor
      $LRGImport::dbCore = $dba;
      
      print STDOUT localtime() . "\tCleaning $lrg_id from " . $dba->dbc()->dbname() . "\n" if ($verbose);
      LRGImport::purge_db($lrg_id,$LRG_COORD_SYSTEM_NAME,$purge);
    }
    # Set the db adaptor in LRGImport back to the core adaptor
    $LRGImport::dbCore = $dbCore;
  }
  
  # Annotate Ensembl genes that overlap this LRG region
  if ($overlap) {
    # Set the db adaptors in the LRGMapping module
    $LRGMapping::dbCore_rw = $dbCore;
    $LRGMapping::dbCore_ro = $dbCore;
    
    #ÊGet a LRG slice
    print STDOUT localtime() . "\tGetting a slice for $lrg_id\n" if ($verbose);
    my $lrg_slice = $sa->fetch_by_region($LRG_COORD_SYSTEM_NAME,$lrg_id) or die("Could not fetch a slice object for " . $LRG_COORD_SYSTEM_NAME . ":" . $lrg_id);
    
    # Get genes that overlap this LRG
    print STDOUT localtime() . "\tGetting genes overlapping $lrg_id\n" if ($verbose);
    my $genes = LRGMapping::get_overlapping_genes($lrg_slice);
    
    # For each overlapping gene, create an XML feature node and check if the overlap is partial or not
    foreach my $gene (@{$genes}) {
      my $feature_node = LRGMapping::gene_2_feature($gene,$lrg_slice);
      print STDOUT localtime() . "\tAdding $lrg_id " . (defined($feature_node->findNode('partial')) ? 'partial ' : '') . "overlap attribute for gene " . $gene->stable_id . " (" . $gene->description . ")\n" if ($verbose);
      LRGImport::add_lrg_overlap($gene->stable_id,$lrg_id,defined($feature_node->findNode('partial')));
    }
  }
  
  if ($import || $add_xrefs || $verify) {
    
  # If lrg_id has been specified but not input_file and a XML file is required, try to fetch it from the LRG website to the /tmp directory
    if (!defined($input_file)) {
    
      print STDOUT localtime() . "\tNo input XML file specified for $lrg_id, attempting to get it from the LRG server\n" if ($verbose);
      my $result = LRGImport::fetch_remote_lrg($lrg_id,[$LRG_EXTERNAL_XML]);
      if ($result->{'success'}) {
	$input_file = $result->{'xmlfile'};
	print STDOUT localtime() . "\tSuccessfully downloaded XML file for $lrg_id and stored it in $input_file\n" if ($verbose);
      }
      else {
	warn("Could not fetch XML file for $lrg_id from external db. Server said: " . $result->{$LRG_EXTERNAL_XML}{'message'});
	warn("Skipping $lrg_id!\n");
	next;
      }
    }
  
    die("ERROR: Input file $input_file does not exist!") unless(-e $input_file);
    
    # create an LRG object from it
    print STDOUT localtime() . "\tCreating LRG object from input XML file $input_file\n" if ($verbose);
    my $lrg = LRG::LRG::newFromFile($input_file) or die("ERROR: Could not create LRG object from XML file!");
    
    # find the LRG ID
    my $lrg_name = $lrg->findNode("fixed_annotation/id")->content();
    print STDOUT localtime() . "\tLRG ID is $lrg_name\n" if ($verbose);
    die("ERROR: Problem with LRG identifier '$lrg_name'") unless ($lrg_name =~ /^LRG\_[0-9]+$/);
    die("ERROR: LRG identifier in $input_file is '$lrg_name' but expected '$lrg_id'") if ($lrg_name ne $lrg_id);
    
    if ($import) {
      
      # Check if the LRG already exists in the database (if the seq_region exists), in which case it should first be deleted
      my $cs_id = LRGImport::add_coord_system($LRG_COORD_SYSTEM_NAME);
      my $seq_region_id = LRGImport::get_seq_region_id($lrg_id,$cs_id);
      if (defined($seq_region_id)) {
	warn("$lrg_id already exists in $coredb\. If you want to replace it, delete it first using the -clean parameter. Will skip it for now");
	#ÊUndefine the input_file so that the next one will be fetched
	undef($input_file);
	#ÊNote, this will also skip xref and verify methods for this LRG
        next;
      }
      
      # Get the assembly that the database uses
      print STDOUT localtime() . "\tGetting assembly name from core db\n" if ($verbose);
      my $db_assembly = LRGImport::get_assembly();
      print STDOUT localtime() . "\tcore db assembly is $db_assembly\n" if ($verbose);
      
      # Find the mapping in the XML file corresponding to the core assembly
      print STDOUT localtime() . "\tGetting mapping from XML file\n" if ($verbose);
      my $updatable_annotation = $lrg->findNode("updatable_annotation");
      my $annotation_sets = $updatable_annotation->findNodeArray("annotation_set");
      my $mapping_node;
      foreach my $annotation_set (@{$annotation_sets}) {
	$mapping_node = $annotation_set->findNode("mapping",{'assembly' => $db_assembly});
	last unless !$mapping_node;
      }
      
      # Warn and skip if the correct mapping could not be fetched
      if (!defined($mapping_node)) {
	warn("Could not find the LRG->Genome mapping corresponding to the core assembly ($db_assembly) for $lrg_id Skipping!");
	#ÊUndefine the input_file so that the next one will be fetched
	undef($input_file);
	#ÊNote, this will also skip xref and verify methods for this LRG
	next;
      }
      
      #ÊWarn if the assembly used is not flagged as the most recent
      warn("Assembly $db_assembly is currently not flagged as the most recent in the $lrg_id XML file!") unless ($mapping_node->{'data'}{'most_recent'} == 1);
      
      my $assembly = $mapping_node->data->{'assembly'};
      print STDOUT localtime() . "\tMapped assembly is $assembly\n" if ($verbose);
      
      # Extract the genomic LRG sequence
      my $lrg_seq = $lrg->findNode('fixed_annotation/sequence')->content();
      #ÊGet the reference genomic sequence from database
      my $chr_name = $mapping_node->data->{'chr_name'};
      my $chr_start = $mapping_node->data->{'start'};
      my $chr_end = $mapping_node->data->{'end'}; 
      my $chr_seq = $sa->fetch_by_region('chromosome',$chr_name,$chr_start,$chr_end)->seq();
      
      # Create pairs array based on the data in the mapping node
      print STDOUT localtime() . "\tCreating pairs from mapping\n" if ($verbose);
      my $mapping = LRGMapping::mapping_2_pairs(
	$mapping_node,
	$lrg_seq,
	$chr_seq
      );
      my $pairs = $mapping->{'pairs'};
      
      # Insert entries for the analysis
      print STDOUT localtime() . "\tAdding analysis data for LRG to $coredb\n" if ($verbose);
      my $analysis_id = LRGImport::add_analysis($LRG_ANALYSIS_LOGIC_NAME);
	
      LRGImport::add_analysis_description(
        $analysis_id,
        $LRG_ANALYSIS_DESCRIPTION,
        $LRG_ANALYSIS_DISPLAY_LABEL,
        0
      );
      
      #ÊLoop over the db adaptors where information will be mirrored and insert in the ones that are defined
      foreach my $dba (@db_adaptors) {
	# Skip this db adaptor if it's not defined
	next unless(defined($dba));
	
	# Set the dbCore variable in LRGImport to the current db adaptor
	$LRGImport::dbCore = $dba;
	
	my $dbname = $dba->dbc()->dbname();
	
	#ÊAdd mapping between the LRG and chromosome coordinate systems to the core db
	print STDOUT localtime() . "\tAdding mapping between $LRG_COORD_SYSTEM_NAME and chromosome coordinate system to $dbname for $lrg_name\n" if ($verbose);
	LRGImport::add_mapping(
	  $lrg_name,
	  $LRG_COORD_SYSTEM_NAME,
	  length($lrg_seq),
	  $mapping
	);
      }
      
      # Set the dbCore variable in LRGImport back to the db adaptor for the core db
      $LRGImport::dbCore = $dbCore;
      
      # Add the transcripts to the core db
      print STDOUT localtime() . "\tAdding transcripts for $lrg_name to core db\n" if ($verbose);
      LRGImport::add_annotation(
	$lrg,
	$lrg_name,
	$LRG_COORD_SYSTEM_NAME,
	$LRG_BIOTYPE,
	$LRG_ANALYSIS_LOGIC_NAME
      );
      
      print STDOUT localtime() . "\tImport done!\n" if ($verbose);
    }
    if ($add_xrefs) {
      
      # This should no longer be done from this script but be included in the main core xref mapping. Will allow it for now.
      #Êdie("Adding xrefs is no longer done from this script. Exiting!");
      
      #ÊGet the Ensembl gene_id for the LRG gene
      my $gene_id = LRGImport::get_object_id_by_stable_id('gene',$lrg_name);
      if (!$gene_id) {
	warn("Could not find gene with stable id $lrg_name in core database! Skipping $lrg_name");
	#ÊUndefine the input_file so that the next one will be fetched
	undef($input_file);
	#ÊNote, this will also skip verify method for this LRG
        next;
      }
      
      # Get the HGNC identifier from the XML 
      my $lrg_gene_name_node = $lrg->findNode("updatable_annotation/annotation_set/lrg_gene_name",{'source' => 'HGNC'});
      if (!$lrg_gene_name_node) {
	warn("Could not find HGNC identifier in XML file! Skipping $lrg_name");
	#ÊUndefine the input_file so that the next one will be fetched
	undef($input_file);
	#ÊNote, this will also skip xref and verify methods for this LRG
        next;
      }
      my $hgnc_name = $lrg_gene_name_node->content();
      
      # A bit cumbersome but.. get the HGNC accession from the XML
      my $annotation_sets = $lrg->findNodeArray('updatable_annotation/annotation_set');
      my $annotation_set_ensembl;
      while (my $annotation_set = shift(@{$annotation_sets})) {
	if ($annotation_set->findNode('source/name')->content() eq 'Ensembl') {
	  $annotation_set_ensembl = $annotation_set;
	  last;
	}
      }
      my $lrg_gene = $annotation_set_ensembl->findNode('features/gene',{'symbol' => $hgnc_name});
      my $hgnc_accession = $lrg_gene->findNode('db_xref',{'source' => 'HGNC'})->data()->{'accession'};
      
      my $xref_id;
      my $object_xref_id;
      
      if (defined($hgnc_accession)) {
	# Add HGNC entry to xref table (or get xref_id if it already exists)
	$xref_id = LRGImport::add_xref('HGNC',$hgnc_accession,$hgnc_name);
	# Add an object_xref for the HGNC xref
	$object_xref_id = LRGImport::add_object_xref($gene_id,'Gene',$xref_id);
      }
      
=head
Skip this, this is done by the external data files 
      #ÊAdd the LRG website as an external db (if not already present)
      LRGImport::add_external_db(
	$LRG_EXTERNAL_DB_NAME,
	$LRG_EXTERNAL_STATUS,
	$LRG_EXTERNAL_PRIORITY,
	$LRG_EXTERNAL_DB_DISPLAY_NAME,
	$LRG_EXTERNAL_DB_RELEASE,
	$LRG_EXTERNAL_DB_ACC_LINKABLE,
	$LRG_EXTERNAL_DB_LABEL_LINKABLE,
	$LRG_EXTERNAL_TYPE
      );
      
      #ÊAdd the Ensembl LRG display as an external db (if not already present). One each for Gene, Transcript
      foreach my $type (('gene','transcript')) {
	LRGImport::add_external_db(
	  $LRG_ENSEMBL_DB_NAME . '_' . $type,
	  $LRG_ENSEMBL_STATUS,
	  $LRG_ENSEMBL_PRIORITY,
	  $LRG_ENSEMBL_DB_DISPLAY_NAME,
	  $LRG_ENSEMBL_DB_RELEASE,
	  $LRG_ENSEMBL_DB_ACC_LINKABLE,
	  $LRG_ENSEMBL_DB_LABEL_LINKABLE,
	  $LRG_ENSEMBL_TYPE
	);
      }
=cut      
      # Add external LRG link to xref table
      my $ext_xref_id = LRGImport::add_xref($LRG_EXTERNAL_DB_NAME,$lrg_name,$lrg_name,undef,'Locus Reference Genomic record for ' . $hgnc_name,'DIRECT');
      
      #ÊAdd an object_xref for the LRG xref
      $object_xref_id = LRGImport::add_object_xref($gene_id,'Gene',$ext_xref_id);
      
      #ÊUpdate the gene table to set the display_xref_id to the LRG xref
      LRGImport::update_rows([qq{display_xref_id = $ext_xref_id}],[qq{gene_id = $gene_id}],['gene']);
      
      #ÊAdd xrefs to the Ensembl coordinate system for the LRG gene
      
      #ÊGet the annotated Ensembl xrefs from the XML file for the LRG gene
      my $lrg_gene_xrefs = $lrg_gene->findNodeArray('db_xref',{'source' => 'Ensembl'});
      
      # Add or get xref_ids for the Ensembl xrefs, the external_db name is Ens_Hs_gene
      foreach my $lrg_gene_xref (@{$lrg_gene_xrefs}) {
	my $stable_id = $lrg_gene_xref->data->{'accession'};
	
	$xref_id = LRGImport::add_xref('Ens_Hs_gene',$stable_id,$stable_id);
	#ÊAdd an object_xref for the LRG xref
	$object_xref_id = LRGImport::add_object_xref($gene_id,'Gene',$xref_id);
	
	my $core_stable_id = $lrg_name;
	
	#ÊDo the same for the Ensembl gene to the Ensembl LRG display
	$xref_id = LRGImport::add_xref($LRG_ENSEMBL_DB_NAME . '_gene',$core_stable_id,$lrg_name);
	# Get the gene_id for the Ensembl gene
	my $core_id = LRGImport::get_object_id_by_stable_id('gene',$stable_id);
	$object_xref_id = LRGImport::add_object_xref($core_id,'Gene',$xref_id);
	# Add xref to the external LRG database
	$object_xref_id = LRGImport::add_object_xref($core_id,'Gene',$ext_xref_id);
      }
      
      # Get Ensembl accessions for transcripts corresponding to transcripts in the fixed section
      my $lrg_transcripts = $lrg_gene->findNodeArray('transcript',{'source' => 'Ensembl'});
      foreach my $lrg_transcript (@{$lrg_transcripts}) {
	my $fixed_id = $lrg_transcript->{'data'}{'fixed_id'};
	my $core_accession = $lrg_transcript->{'data'}{'transcript_id'};
	next unless(defined($fixed_id) && defined($core_accession));
	
	# Get the core db LRG transcript_id for this transcript
	my $core_stable_id = $lrg_name . '_' . $fixed_id;
	my $core_id = LRGImport::get_object_id_by_stable_id('transcript',$core_stable_id);
	next unless(defined($core_id));
	
	$xref_id = LRGImport::add_xref('Ens_Hs_transcript',$core_accession,$core_accession);
	#ÊAdd an object_xref for the LRG xref
	$object_xref_id = LRGImport::add_object_xref($core_id,'Transcript',$xref_id);
	
	#ÊDo the same for the Ensembl transcript to the Ensembl LRG display
	$xref_id = LRGImport::add_xref($LRG_ENSEMBL_DB_NAME . '_transcript',$core_stable_id,$core_stable_id);
	# Get the gene_id for the Ensembl gene
	my $core_id = LRGImport::get_object_id_by_stable_id('transcript',$core_accession);
	$object_xref_id = LRGImport::add_object_xref($core_id,'Transcript',$xref_id);
	
	# Do the same for the translation
	my $lrg_protein = $lrg_transcript->findNode('protein_product',{'source' => 'Ensembl'});
	next unless(defined($lrg_protein));
	$core_accession = $lrg_protein->{'data'}{'accession'};
	next unless(defined($core_accession));
	$core_id = LRGImport::get_translation_id($core_id);
	
	$xref_id = LRGImport::add_xref('Ens_Hs_translation',$core_accession,$core_accession);
	#ÊAdd an object_xref for the LRG xref
	$object_xref_id = LRGImport::add_object_xref($core_id,'Translation',$xref_id);
	
      }
      
    }
    
    #ÊCheck that the mapping stored in the database give the same sequences as those stored in the XML file
    if ($verify) {
      
      # Needs to flush the db adaptor when doing an import before verifying content. Not sure how to do this, so will only print a warning
      warn("Verifying the import now will probably tell you that there are inconsistencies because the db adaptor haven't been flushed. You should re-run the script doing verification after the import has finished") if ($import);
      
      #ÊA flag to inddicate if everything is ok
      my $passed = 1;
      
      # Get the genomic sequence from the XML file
      my $genomic_seq_xml = $lrg->findNode("fixed_annotation/sequence")->content();
      # Get a slice from the database corresponding to the LRG
      my $lrg_slice = $sa->fetch_by_region($LRG_COORD_SYSTEM_NAME,$lrg_id);
      if (!defined($lrg_slice)) {
	warn("Could not fetch a slice object for " . $LRG_COORD_SYSTEM_NAME . ":" . $lrg_id);
	$passed = 0;
      }
      else {
	# Get the genomic sequence of the slice
	my $genomic_seq_db = $lrg_slice->seq();
	
	# Compare the sequences
	if ($genomic_seq_xml ne $genomic_seq_db) {
	  warn("Genomic sequence from core db is different from genomic sequence in XML file for $lrg_id");
	  $passed = 0;  
	}
	
	#ÊCompare each transcript
	my $transcripts_xml = $lrg->findNodeArray('fixed_annotation/transcript');
	my $transcripts_db = $lrg_slice->get_all_Transcripts(undef,$LRG_ANALYSIS_LOGIC_NAME);
	foreach my $transcript_xml (@{$transcripts_xml}) {
	  # Get the fixed id
	  my $fixed_id = $transcript_xml->{'data'}{'name'};
	  #ÊThe expected transcript_stable_id based on the XML fixed id
	  my $stable_id = $lrg_id . '_' . $fixed_id;
	  # Get the ensembl transcript with the corresponding stable_id
	  my @db_tr = grep {$_->stable_id() eq $stable_id} @{$transcripts_db};
	  # Check that we got one transcript back
	  if (!defined(@db_tr) || scalar(@db_tr) != 1) {
	    warn("Could not unambiguously get the correct Ensembl transcript corresponding to $lrg_id $fixed_id");
	    $passed = 0;
	    next;
	  }
	  
	  my $transcript_db = $db_tr[0];
	  
	  #ÊGet the cDNA sequence from the XML file
	  my $cDNA_xml = $transcript_xml->findNode('cdna/sequence')->content();
	  # Get the cDNA sequence from the db
	  my $cDNA_db = $transcript_db->spliced_seq();
	  # Compare the sequences
	  if ($cDNA_xml ne $cDNA_db) {
	    warn("cDNA sequence from core db is different from cDNA sequence in XML file for $lrg_id transcript $fixed_id");
	    $passed = 0;
	    next;
	  }
	  
	  #ÊGet the translation from the XML file
	  my $translation_xml = $transcript_xml->findNode('coding_region/translation/sequence')->content();
	  #ÊGet the translation from the db
	  my $translation_db = $transcript_db->translation()->seq();
	  # Remove any terminal stop codons
	  $translation_xml =~ s/\*$//;
	  $translation_db =~ s/\*$//;
	  
	  # Compare the sequences
	  if ($translation_xml ne $translation_db) {
	    warn("Peptide sequence from core db is different from peptide sequence in XML file for $lrg_id transcript $fixed_id");
	    $passed = 0;
	    next;
	  }	
	}
      }
      
      if ($passed) {
	print STDOUT "$lrg_id is consistent between XML file and core db\n";
      }
      else {
	print STDOUT "$lrg_id has inconsistencies between XML file and core db\n";
      }
    }
  }
  #ÊUndefine the input_file so that the next one will be fetched
  undef($input_file);
}

sub usage {
	
  print qq{
  Usage: perl import.lrg.pl [OPTION]
  
  Import or update/remove a LRG record in a Core database
	
  Options:
    
    Database credentials are specified on the command line
    
      -host		Core database host name (Required)
      -port		Core database port (Required)
      -dbname		Core database name (Required)
      -user		Core database user (Required)
      -pass		Core database password (Optional)
    
    In addition, names for otherfeatures, cdna and vega databases can be specified on the command line. If so,
    the mapping data will be inserted into those databases as well. They are assumed to reside on the same host as
    the core database and be accessed with the same credentials. Before doing an import, will try to sync the
    analysis, coord_system, meta and seq_region tables w.r.t. their Auto_increment values.
    
      -otherfeatures	Otherfeatures database name
      -cdna		cDNA database name
      -vega		vega database name
      
    An input file can be specified. This is required when reverting the Core database. If an input file is
    specified when importing, verifying, cleaning, adding xrefs or annotating overlaps, all specified LRG
    identifiers are overridden and only the LRG in the input XML file is processed.
    
      -input_file	LRG XML file when importing or adding xrefs
			Tab-separated file with table, field and max-values when reverting database
			to a previous state.
			
    Any number of LRG identifiers can be specified. Each LRG will then be processed in turn. If an identifier is
    specified when importing, verifying or adding xrefs, the script will attempt to download the corresponding XML
    file from the LRG website.
    
      -lrg_id		LRG identifier on the form LRG_N, where N is an integer
      
    If neither input file nor lrg identifiers are specified, the script will obtain a list of publicly available
    LRGs from the LRG ftp site and the user can interactively choose which LRGs to process.
    
    What action the script will perform is dictated by the following flags:
    
      -import		The data in the supplied, or downloaded, LRG XML file will be imported into the Core
			database. This includes the mapping to the current assembly and the
			transcripts in the fixed section
      
      -verify		Verify the consistency between the sequences stored in the LRG XML file and what the
			API gets when accessing the core db. Will check genomic sequence, cDNA (spliced transcript)
			and peptide translation.
      
      -xrefs		This will add xrefs to HGNC, the external LRG website as well as to the
			corresponding Ensembl genes for genes on the lrg coordinate system.
			Will also add xrefs to the corresponding Ensembl genes, linking them to
			the external LRG website and to the genes on the lrg coordinate system.
			This should only need to be run for release 58. Subsequent releases will
			take care of this through the normal xref pipeline, using HGNC data
			
      -overlap		For the LRG specified by the -lrg_id argument, find all Ensembl genes in the
			chromosome coordinate system that overlap. Add a gene_attrib to these, indicating
			that they overlap a LRG region, either partially or completely. Note that the LRG
			must already be stored in the Core database
      
      -clean		Remove all entries in the Core database specifically relating to the
			LRG that was specified with the -lrg_id argument
			
      -purge		If specified, will remove EVERYTHING LRG related from the database(s) (e.g. coord_system, analysis)
			once all LRG seq_regions have been removed. Will not remove anything if seq_regions still exist.
			
      -max		Dump a tab-separated list of table, field and max-values for tables
			affected by a LRG import to stdout. If run before an import, this data can
			be used to revert the database after an import
			
      -revert		Don't use this unless you are really sure of what it does! Will delete all rows having
			field values greater than the max-value supplied via the tab-separated input file. These
			should have been generated using the -max mode. Beware that this will delete all entries
			added after the -max command was run, not necessarily just your own!
			
      -verbose		Progress information is printed
      -help		Print this message
      
  };
  exit(0);
}

