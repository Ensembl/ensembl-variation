#!perl

use strict;

use LRG;

# get the files to read from and write to as three command-line arguments
my ($file1, $file2, $outfile, $merge_mapping) = @ARGV;

# sanity checks
die "ERROR: Could not read from input file $file1\n" unless -e $file1;
die "ERROR: Could not read from input file $file2\n" unless -e $file2;
die "ERROR: No output file specified\n" unless defined $outfile;

# create LRG objects for each of the input files
my $root1 = LRG::LRG::newFromFile($file1, $outfile);
my $root2 = LRG::LRG::newFromFile($file2);

# Sanity check that each has the same LRG ID
# The findNode() method searches the XML tree for nodes - separating by "/"s represents parent/child heirarchy
# such that in our example id should be a child node (although not necessary an immediate child) of the
# fixed_annotation node
# We use the content() method to get the plain-text content of a node
my $id1 = $root1->findNode("fixed_annotation/id")->content();
my $id2 = $root2->findNode("fixed_annotation/id")->content();

die "ERROR: LRG identifiers do not match - file 1 is $id1, file 2 is $id2\n" unless $id1 eq $id2;

my $set2 = $root2->findNode("annotation_set");

die "ERROR: Could not find annotation set in $file2\n" unless defined $set2;

my $insert_point = $root1->findNode("updatable_annotation");

die "ERROR: Could not find updatable_annotation node in $file1\n" unless defined $insert_point;

# we can use the addExisting method to add an existing node and associated child nodes to another node tree
$insert_point->addExisting($set2);

# This is to create a separate annotation_set for shared mapping data and gene symbol
if ($merge_mapping) {
    my $newnode = LRG::LRG::new();
    $newnode = $newnode->addNode("annotation_set");
    my $srcnode = $newnode->addNode("source");
    $srcnode->addNode("name")->content("LRG");
    $srcnode->addNode("url")->content("http://www.lrg-sequence.org/");
    $srcnode = $srcnode->addNode("contact");
    $srcnode->addNode("name")->content("Locus Reference Genomic");
    $srcnode->addNode("url")->content("http://www.lrg-sequence.org/page.php?page=contact");
    $srcnode->addNode("email")->content("feedback\@lrg-sequence.org");
    $newnode->addNode("modification_date")->content($root1->date);
    
 # Get the mapping details
    $srcnode = $root1->findNodeArray("updatable_annotation/annotation_set/mapping");
    if ($srcnode && scalar(@{$srcnode}) > 0) {
        foreach my $src (@{$srcnode}) { 
            $newnode->addExisting($src);
        }
    }
    $srcnode = $root2->findNodeArray("updatable_annotation/annotation_set/mapping");
    if ($srcnode && scalar(@{$srcnode}) > 0) {
        foreach my $src (@{$srcnode}) { 
            $newnode->addExisting($src);
        }
    }
    $srcnode = $root1->findNode("updatable_annotation/annotation_set/lrg_gene_name");
    $newnode->addExisting($srcnode) unless (!$srcnode);
    $srcnode = $root2->findNode("updatable_annotation/annotation_set/lrg_gene_name");
    $newnode->addExisting($srcnode) unless (!$srcnode);
    
    $insert_point->addExisting($newnode);
}

# the printAll method then writes the structure to XML; the filename for the XML
# is specified as the second argument in the new() method call above. If no output
# file is given in this argument, a temporary file is used instead.
$root1->printAll();
