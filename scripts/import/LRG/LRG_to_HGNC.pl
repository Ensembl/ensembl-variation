#!perl

use strict;

use Getopt::Long;
use LRG::LRG;


# Default options
my $lrg_directory;
my $lrg_name;
my $lrg_file;
my $help;

usage() unless scalar @ARGV;

# Get options from command line
GetOptions(
    'lrg_directory=s' => \$lrg_directory,
    'lrg_name=s' => \$lrg_name,
    'lrg_file=s' => \$lrg_file,
    'help!' => \$help
);

usage() if($help);

die("Either the name of an LRG or the path to the LRG XML file must be specified!") unless (defined($lrg_name) || defined($lrg_file));

$lrg_file ||= $lrg_directory . "/" . $lrg_name . ".xml";

die("File $lrg_file does not exist!") unless (-e $lrg_file);

#ÊCreate an LRG object from the input file
my $root = LRG::LRG::newFromFile($lrg_file);

# Extract the LRG ID
my $lrg_id = $root->findNode("fixed_annotation/id")->content();

# Raise a warning in case the identifier is not in the expected format
warn("LRG identifier '$lrg_id' is not in the expected format (LRG_[0-9]+)") unless ($lrg_id =~ /^LRG\_[0-9]+$/);

# Raise a warning if the LRG identifier within the XML file is different from what was specified on the command line
warn("LRG identifier '$lrg_id', found in the XML file '$lrg_file', is different from what was specified on the command line ('$lrg_name')") unless (!defined($lrg_name) || $lrg_name eq $lrg_id);

#ÊExtract the annotation_sets from the updatable annotation section
my $annotation_sets = $root->findNodeArray('updatable_annotation/annotation_set');

my $hgnc_symbol;

# Loop over the annotation sets and get the HGNC gene symbol
while (my $annotation_set = shift @{$annotation_sets}) {
    my $source_node = $annotation_set->findNode('source/name');
    next unless (defined($source_node) && $source_node->content() eq 'LRG');
    
    my $hgnc_node = $annotation_set->findNode('lrg_gene_name',{'source' => 'HGNC'});
    $hgnc_symbol = $hgnc_node->content() unless (!defined($hgnc_node));
}

die("Could not find HGNC gene symbol within LRG section of XML file!") unless (defined($hgnc_symbol));

#ÊOutput the LRG to HGNC mapping
print STDOUT $lrg_id . "\t" . $hgnc_symbol . "\n";


#ÊDisplay usage and then exit
sub usage() {
	
    print STDOUT qq{
    Usage:
    
      perl LRG_to_HGNC.pl [options] INPUT 
    
      Get the HGNC symbol corresponding to a LRG record
    
      INPUT: Either -lrg_name or -lrg_file must be specified
        -lrg_name           LRG id
        -lrg_file           Path to XML file
    
      Options:
        -lrg_directory	path to directory containing XML file corresponding to specified LRG
        -h, --help		print this message
    };
    
    die;
}


