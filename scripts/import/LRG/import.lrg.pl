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

#!/software/bin/perl

use strict;

use Getopt::Long;
use List::Util qw (min max);
use LWP::UserAgent;
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
my $LRG_ANALYSIS_WEB_DATA = qq{{'colour_key' => 'rna_[status]','caption' => 'LRG gene','label_key' => '[text_label] [display_label]','name' => 'LRG Genes','default' => {'MultiTop' => 'gene_label','contigviewbottom' => 'transcript_label','MultiBottom' => 'collapsed_label','contigviewtop' => 'gene_label','alignsliceviewbottom' => 'as_collapsed_label','cytoview' => 'gene_label'},'multi_caption' => 'LRG genes','key' => 'ensembl'}};
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
my $LRG_EXTERNAL_XML = q{ftp://ftp.ebi.ac.uk/pub/databases/lrgex/};#http://www.lrg-sequence.org/lrg.php};

my $host;
my $port;
my $dbname;
my $user;
my $pass;
my $help;
my $clean;
my $verbose;
my $overlap;
my $lrg_id;
my $input_file;
my $import;
my $add_xrefs;
my $max_values;
my $revert;

usage() if (!scalar(@ARGV));

# get options from command line
GetOptions(
  'host=s'		=> \$host,
  'port=i'		=> \$port,
  'dbname=s'		=> \$dbname,
  'user=s'		=> \$user,
  'pass=s'		=> \$pass,
  'help!' 		=> \$help,
  'verbose!' 		=> \$verbose,
  'clean!' 		=> \$clean,
  'overlap!' 		=> \$overlap,
  'lrg_id=s' 		=> \$lrg_id,
  'input_file=s' 	=> \$input_file,
  'import!' 		=> \$import,
  'xrefs!' 		=> \$add_xrefs,
  'max!' 		=> \$max_values,
  'revert!' 		=> \$revert
);

usage() if (defined($help));

die("Database credentials (-host, -port, -dbname, -user) need to be specified!") unless (defined($host) && defined($port) && defined($dbname) && defined($user));
die("Supplied LRG id is not in the correct format ('LRG_NNN')") if (defined($lrg_id) && $lrg_id !~ m/^LRG\_[0-9]+$/);
die("A LRG id needs to be specified via -lrg_id parameter!") if (!defined($lrg_id) && ($clean || $overlap));

# If lrg_id has been specified but not input_file and a XML file is required, try to fetch it from the LRG website to the /tmp directory
if ((defined($import) || defined($add_xrefs)) && !defined($input_file) && defined($lrg_id)) {
  
  print STDOUT localtime() . "\tNo input XML file specified for $lrg_id, attempting to get it from the LRG server\n" if ($verbose);
  
  # Create a user agent object
  my $ua = LWP::UserAgent->new;

  # Create a request
  my $req = HTTP::Request->new(GET => "$LRG_EXTERNAL_XML/$lrg_id\.xml");
  #$req->content("id=$lrg_id");

  # Pass request to the user agent and get a response back
  my $res = $ua->request($req);

  # Check the outcome of the response
  die("Could not fetch XML file for $lrg_id from server. Server says: " . $res->message) unless ($res->is_success);
  
  # Write the XML to a temporary file
  $input_file = '/tmp/' . $lrg_id . '.xml';
  open(XML_OUT,'>',$input_file) or die("Could not open $input_file for writing");
  print XML_OUT $res->content;
  close(XML_OUT);
  
  print STDOUT localtime() . "\tSuccessfully downloaded XML file for $lrg_id and stored it in $input_file\n" if ($verbose);
}

die("An input LRG XML file (or a LRG identifier to fetch it) must be specified in order to import into core db!") if ((defined($import) || defined($add_xrefs)) && !defined($input_file));
die("A tab-separated input file with table, field and max_value columns must be specified in order to revert the core db!") if (defined($revert) && !defined($input_file));

# Connect to core database
print STDOUT localtime() . "\tGetting human core db adaptor\n" if ($verbose);
my $dbCore = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -host => $host,
  -user => $user,
  -pass => $pass,
  -port => $port,
  -dbname => $dbname
) or die("Could not get a database adaptor to $dbname on $host:$port");
print STDOUT localtime() . "\tConnected to $dbname on $host:$port\n" if ($verbose);

$LRGImport::dbCore = $dbCore;

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

# Clean up data in the database if required
if ($clean) {
  print STDOUT localtime() . "\tCleaning $clean from core db\n" if ($verbose);
  LRGImport::purge_db($lrg_id,$LRG_COORD_SYSTEM_NAME);
}
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
    print STDOUT localtime() . "\tAdding $lrg_id " . (defined($feature_node->findNode('partial')) ? 'partial ' : '') . "overlap attribute for gene $gene->stable_id ($gene->description)\n" if ($verbose);
    LRGImport::add_lrg_overlap($gene->stable_id,$lrg_id,defined($feature_node->findNode('partial')));
  }
}
if ($import || $add_xrefs) {
  die("ERROR: Input file $input_file does not exist!") unless(-e $input_file);
  
  # create an LRG object from it
  print STDOUT localtime() . "\tCreating LRG object from input XML file $input_file\n" if ($verbose);
  my $lrg = LRG::LRG::newFromFile($input_file) or die("ERROR: Could not create LRG object from XML file!");
  
  # find the LRG ID
  my $lrg_name = $lrg->findNode("fixed_annotation/id")->content();
  print STDOUT localtime() . "\tLRG ID is $lrg_name\n" if ($verbose);
  die("ERROR: Problem with LRG identifier '$lrg_name'") unless ($lrg_name =~ /^LRG\_[0-9]+$/);
  
  if ($import) {
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
    
    # Die if the correct mapping could not be fetched
    die("Could not find the LRG->Genome mapping corresponding to the core assembly ($db_assembly)") unless (defined($mapping_node));
    
    #ÊWarn if the assembly used is not flagged as the most recent
    warn("Assembly $db_assembly is currently not flagged as the most recent in the XML file!") unless ($mapping_node->{'data'}{'most_recent'} == 1);
    
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
    print STDOUT localtime() . "\tAdding analysis data for LRG to core db\n" if ($verbose);
    my $analysis_id = LRGImport::add_analysis($LRG_ANALYSIS_LOGIC_NAME);
    
    LRGImport::add_analysis_description(
      $analysis_id,
      $LRG_ANALYSIS_DESCRIPTION,
      $LRG_ANALYSIS_DISPLAY_LABEL,
      1,
      $LRG_ANALYSIS_WEB_DATA
    );
    
    #ÊAdd mapping between the LRG and chromosome coordinate systems to the core db
    print STDOUT localtime() . "\tAdding mapping between $LRG_COORD_SYSTEM_NAME and chromosome coordinate system to core db for $lrg_name\n" if ($verbose);
    LRGImport::add_mapping(
      $lrg_name,
      $LRG_COORD_SYSTEM_NAME,
      length($lrg_seq),
      $mapping
    );
    
    # Add the transcripts to the core db
    print STDOUT localtime() . "\tAdding transcripts for $lrg_name to core db\n" if ($verbose);
    LRGImport::add_annotation(
      $lrg,
      $lrg_name,
      $LRG_COORD_SYSTEM_NAME,
      $LRG_BIOTYPE,
      $LRG_ANALYSIS_LOGIC_NAME
    );
    
    print STDOUT localtime() . "\tAll done, exiting!\n" if ($verbose);
  }
  if ($add_xrefs) {
    
    #ÊGet the Ensembl gene_id for the LRG gene
    my $gene_id = LRGImport::get_object_id_by_stable_id('gene',$lrg_name . '_g1') or die ("Could not find gene with stable id $lrg_name\_g1 in core database!");
    
    # Get the HGNC identifier from the XML 
    my $lrg_gene_name_node = $lrg->findNode("updatable_annotation/annotation_set/lrg_gene_name",{'source' => 'HGNC'}) or die ("Could not find HGNC identifier in XML file!");
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
    
    # Add HGNC entry to xref table (or get xref_id if it already exists)
    my $xref_id = LRGImport::add_xref('HGNC',$hgnc_accession,$hgnc_name);
    
    # Add an object_xref for the HGNC xref
    my $object_xref_id = LRGImport::add_object_xref($gene_id,'Gene',$xref_id);
    
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
    
    # Add external LRG link to xref table
    $xref_id = LRGImport::add_xref($LRG_EXTERNAL_DB_NAME,$lrg_name,$lrg_name,undef,'Locus Reference Genomic record for ' . $hgnc_name,'DIRECT');
    
    #ÊAdd an object_xref for the LRG xref
    $object_xref_id = LRGImport::add_object_xref($gene_id,'Gene',$xref_id);
    
    #ÊUpdate the gene table to set the display_xref_id to the LRG xref
    LRGImport::update_rows([qq{display_xref_id = $xref_id}],[qq{gene_id = $gene_id}],['gene']);
    
    #ÊAdd xrefs to the Ensembl coordinate system for the LRG gene
    
    #ÊGet the annotated Ensembl xrefs from the XML file for the LRG gene
    my $lrg_gene_xrefs = $lrg_gene->findNodeArray('db_xref',{'source' => 'Ensembl'});
    
    # Add or get xref_ids for the Ensembl xrefs, the external_db name is Ens_Hs_gene
    foreach my $lrg_gene_xref (@{$lrg_gene_xrefs}) {
      my $stable_id = $lrg_gene_xref->data->{'accession'};
      $xref_id = LRGImport::add_xref('Ens_Hs_gene',$stable_id,$stable_id);
      #ÊAdd an object_xref for the LRG xref
      $object_xref_id = LRGImport::add_object_xref($gene_id,'Gene',$xref_id);
      
      my $core_stable_id = $lrg_name . '_g1';
      
      #ÊDo the same for the Ensembl gene to the Ensembl LRG display
      $xref_id = LRGImport::add_xref($LRG_ENSEMBL_DB_NAME . '_gene',$core_stable_id,$lrg_name);
      # Get the gene_id for the Ensembl gene
      my $core_id = LRGImport::get_object_id_by_stable_id('gene',$stable_id);
      $object_xref_id = LRGImport::add_object_xref($core_id,'Gene',$xref_id);
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
}

sub usage() {
	
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
      
    A LRG identifier must be specified when cleaning a LRG from the database or when annotating
    Ensembl genes that overlap LRG regions. If an identifier is specified when importing or adding
    xrefs, the script will attempt to download the corresponding XML file from the LRG website.
    
      -lrg_id		LRG identifier on the form LRG_N, where N is an integer
      
    An input file can be specified. This is required when importing a LRG (unless the XML file can
    be downloaded from the LRG website, using a supplied identifier) or when reverting the Core database
    
      -input_file	LRG XML file when importing or adding xrefs
			Tab-separated file with table, field and max-values when reverting database
			to a previous state.
			
      
    What action the script will perform is dictated by the following flags:
    
      -import		The data in the supplied, or downloaded, LRG XML file will be imported into the Core
			database. This includes the mapping to the current assembly and the
			transcripts in the fixed section
			
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

