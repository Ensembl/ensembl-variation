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

# set default options
my $registry_file = "ensembl.registry";
my ($help);
my $clean;
my $verbose;
my $overlap;

usage() unless scalar @ARGV;

# get options from command line
GetOptions(
	   'registry_file=s' => \$registry_file,
	   'help' => \$help,
	   'verbose!' => \$verbose,
	   'clean=s' => \$clean,
	   'overlap=s' => \$overlap
);

usage() if $help;

# Connect to core database
print STDOUT localtime() . "\tLoading registry from $registry_file\n" if ($verbose);
Bio::EnsEMBL::Registry->load_all( $registry_file );

print STDOUT localtime() . "\tGetting human core db adaptor\n" if ($verbose);
my $dbCore = Bio::EnsEMBL::Registry->get_DBAdaptor('human','core_rw') or die("ERROR: Could not connect to core database");

$LRGImport::dbCore = $dbCore;

# Clean up data in the database if required
if ($clean) {
  print STDOUT localtime() . "\tCleaning $clean from core db\n" if ($verbose);
  LRGImport::purge_db($clean,$LRG_COORD_SYSTEM_NAME);
  exit(0);
}
if ($overlap) {
  my $gene_attrib = LRGMapping::get_overlapping_genes($overlap,$overlap,$dbCore);
  while (my ($id,$attrib) = each(%{$gene_attrib})) {
    print $id . "\t" . $attrib . "\n";
  }
  exit(0);
}
# check for an input file
my $input_file = shift @ARGV;
die("ERROR: No input file specified or file does not exist!") unless(-e $input_file);

# create an LRG object from it
print STDOUT localtime() . "\tCreating LRG object from input XML file $input_file\n" if ($verbose);
my $lrg = LRG::LRG::newFromFile($input_file) or die("ERROR: Could not create LRG object from XML file!");

# find the LRG ID
my $lrg_name = $lrg->findNode("fixed_annotation/id")->content();
print STDOUT localtime() . "\tLRG ID is $lrg_name\n" if ($verbose);
die("ERROR: Problem with LRG identifier '$lrg_name'") unless ($lrg_name =~ /^LRG\_[0-9]+$/);

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
my $sa = $dbCore->get_SliceAdaptor();
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
  $assembly,
  $chr_name,
  $pairs
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

sub clean {
  my $lrg_name = shift;
  my $lrg_coord_system_name = shift;
  my $dbCore = shift;
  
  $lrg_coord_system_name = $lrg_name;
  
  # Clean out entries from the meta table
  my $stmt = qq{
    DELETE FROM
      meta
    WHERE
      meta_key = 'assembly.mapping' AND
      meta_value LIKE '$lrg_coord_system_name\#\%'
  };
  $dbCore->dbc->do($stmt);
  
  my $cs_id = LRGImport::get_coord_system_id($lrg_coord_system_name,1);
  my $seq_region_id = LRGImport::get_seq_region_id($lrg_name,$cs_id);

  # Remove entry from seq_region    
  $stmt = qq{
    DELETE FROM
      seq_region
    WHERE
      seq_region_id = '$seq_region_id'
  };
  $dbCore->dbc->do($stmt);
    
  # Remove entries from seq_region_attrib table
  $stmt = qq{
    DELETE FROM
      seq_region_attrib
    WHERE
      seq_region_id = '$seq_region_id'
  };
  $dbCore->dbc->do($stmt);
  
  # Remove mapping from assembly table
  $stmt = qq{
    DELETE FROM
      assembly
    WHERE
      asm_seq_region_id = '$seq_region_id'
  };
  $dbCore->dbc->do($stmt);
    
  my $transcript_ids;
  
  # Remove entries from exon, gene and transcript tables that lie on the seq region corresponding to this LRG entry   
  foreach my $table (('exon','gene','transcript')) {
    my $key = $table . '_id';
    my $id_table = $table . '_stable_id';
    
    # Get all concerned ids  
    $stmt = qq{
      SELECT
	$key
      FROM
	$table
      WHERE
	seq_region_id = '$seq_region_id'
    };
    my $ids = $dbCore->dbc->db_handle->selectall_arrayref($stmt);
    
    # Save transcript ids for later as we need them to clear away translations  
    if ($table eq 'transcript') {
      $transcript_ids = \@{$ids};
    }
      
    foreach my $id (@{$ids}) {
      # Remove the stable_id entry
      $stmt = qq{
        DELETE FROM
          $id_table
        WHERE
          $key = '$id->[0]'
      };
      $dbCore->dbc->do($stmt);
      
      # Remove from table
      $stmt = qq{
        DELETE FROM
          $table
        WHERE
          $key = '$id->[0]'
      };
      $dbCore->dbc->do($stmt);
    }
  }
    
  # Drop entries from the translation table
  foreach my $id (@{$transcript_ids}) {
    $stmt = qq{
      SELECT
        translation_id
      FROM
        translation
      WHERE
        transcript_id = '$id->[0]'
      LIMIT 1
    };
    my $translation_id = $dbCore->dbc->db_handle->selectall_arrayref($stmt)->[0][0];
    
    $stmt = qq{
	 DELETE FROM
	   translation_stable_id
	 WHERE
	   translation_id = '$translation_id'
    };
    $dbCore->dbc->do($stmt);
    $stmt = qq{
	 DELETE FROM
	   translation
	 WHERE
	   translation_id = '$translation_id'
    };
    $dbCore->dbc->do($stmt);      
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

