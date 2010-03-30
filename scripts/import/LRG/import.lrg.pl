#!/software/bin/perl

use strict;

use Getopt::Long;
use List::Util qw (min max);
use LRG::LRG;
use LRG::LRGImport;
use LRG::LRGMapping;

use Bio::Seq;
use Bio::EnsEMBL::Registry;

# Some constants
my $LRG_COORD_SYSTEM_NAME = q{LRG};
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
my $LRG_WEBSITE_XREF_ROOT_URL = q{ftp://ftp://ftp.ebi.ac.uk/pub/databases/lrgex/};


my $registry_file;
my $help;
my $clean;
my $verbose;
my $overlap;
my $lrg_id;
my $input_file;
my $import;
my $add_xrefs;

usage() unless scalar @ARGV;

# get options from command line
GetOptions(
  'registry_file=s' => \$registry_file,
  'help' => \$help,
  'verbose!' => \$verbose,
  'clean!' => \$clean,
  'overlap!' => \$overlap,
  'lrg_id=s' => \$lrg_id,
  'input_file=s' => \$input_file,
  'import!' => \$import,
  'xrefs!' => \$add_xrefs
);

usage() if $help;

die("Supplied LRG id is not in the correct format ('LRG_NNN')") if (defined($lrg_id) && $lrg_id !~ m/^LRG\_[0-9]+$/);
die("A LRG id needs to be specified via -lrg_id parameter!") if (!defined($lrg_id) && ($clean || $overlap));
die("An input LRG XML file must be specified in order to import into core db!") if (defined($import) && !defined($input_file));

# Connect to core database
print STDOUT localtime() . "\tLoading registry from $registry_file\n" if ($verbose);
Bio::EnsEMBL::Registry->load_all( $registry_file );

print STDOUT localtime() . "\tGetting human core db adaptor\n" if ($verbose);
my $dbCore = Bio::EnsEMBL::Registry->get_DBAdaptor('human','core_rw') or die("ERROR: Could not connect to core database");
$LRGImport::dbCore = $dbCore;

#ÊGet a slice adaptor
print STDOUT localtime() . "\tGetting slice adaptor\n" if ($verbose);
my $sa = $dbCore->get_SliceAdaptor();

# Clean up data in the database if required
if ($clean) {
  print STDOUT localtime() . "\tCleaning $clean from core db\n" if ($verbose);
  LRGImport::purge_db($lrg_id,$LRG_COORD_SYSTEM_NAME);
}
if ($overlap) {
  # Set the db adaptors in the LRGMapping module
  $LRGMapping::dbCore_rw = $dbCore;
  
  print STDOUT localtime() . "\tGetting human read-only core db adaptor\n" if ($verbose);
  $LRGMapping::dbCore_ro = Bio::EnsEMBL::Registry->get_DBAdaptor('human','core_ro') or die("ERROR: Could not connect to read-only core database");
  
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
    
    # Add to xref table
    $xref_id = LRGImport::add_xref($LRG_EXTERNAL_DB_NAME,$lrg_name . '.xml',$lrg_name,undef,'Locus Reference Genomic record for ' . $hgnc_name,'DIRECT');
    
    #ÊAdd an object_xref for the LRG xref
    $object_xref_id = LRGImport::add_object_xref($gene_id,'Gene',$xref_id);
    
    #ÊUpdate the gene table to set the display_xref_id to the LRG xref
    LRGImport::update_rows({'display_xref_id' => $xref_id},{'gene_id' => $gene_id},['gene']);
  }
}

sub usage() {
	
	print
	"Usage:
	
	perl import.lrg.pl [options] LRG.xml
	
	-r || --registry_file	registry file
	-c || --clean LRG name	remove entries for this LRG from core db
	-v || --verbose		display progress information
	-h || --help		print this message\n";
	
	die;
}

