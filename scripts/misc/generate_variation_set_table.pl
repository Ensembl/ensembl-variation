#! perl -w

#ÊScript to dump out a table of variation sets that can be used in the documentation

use strict;
use warnings;
use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

# Print the usage instructions if run without parameters
usage() unless (scalar(@ARGV));

my $species    = shift;
my $host       = shift;
my $db_version = shift;

die ("Species, db_host and db_version must be specified") unless ($species && $host && $db_version);

# Filters
my @filters = ('fail_');

# Load the registry from db
$registry->load_registry_from_db(
    -host => $host,
    -user => 'ensro',
		-db_version => $db_version
);

# Get a VariationSetAdaptor on the human variation database
my $vs_adaptor = $registry->get_adaptor($species,'variation','variationset');

# Get all top-level variation sets
my $top_vss = $vs_adaptor->fetch_all_top_VariationSets();


my $table_header = qq{
	<tr>
		<th>Name</th>
		<th>Short name</th>
		<th>Description</th>
	</tr>
};
	
# Loop over the top level variation sets and recursively print the subsets
my $com_rowcount = 0;
my $rowcount     = 0;
my $com_sets;
my $sets;
foreach my $top_vs (@{$top_vss}) {
	my $is_com = 0;
	# Common set
	foreach my $com_filter (@filters) {
		if ($top_vs->short_name =~ /^$com_filter/) {
			$com_sets->{$top_vs->short_name} = $top_vs;
			$is_com = 1;
			last;
		}
	}
	# Human specific set
	if (!$is_com) {
		$sets->{$top_vs->short_name} = $top_vs;
	}
}


## Print the common table headers
print "<h4>Variation sets common to all species</h4>\n";
print "<table id=\"variation_set_table\" class=\"ss\">\n";
print "$table_header\n";

foreach my $com_set_name (sort(keys(%$com_sets))) {
	print_set($com_sets->{$com_set_name},\$com_rowcount);
}
print "</table>\n";


## Print the human specific table headers
print "<br />\n<h4>Variation sets specific to Human</h4>\n";
print "<table id=\"human_variation_set_table\" class=\"ss\">\n";
print $table_header;

foreach my $set_name (sort(keys(%$sets))) {
	print_set($sets->{$set_name},\$rowcount);
}
print "</table>\n";



# We define a function that will help us recurse over the set hierarchy and print the data   
sub print_set {
	my $set = shift;
	my $rowcount = shift;
	my $indent = shift || 0;
	
	# Highlight even row numbers
	${$rowcount}++;
	my $rowclass = (${$rowcount}%2 == 0 ? " class=\"bg2\"" : "");
	
	# Put a bullet next to subsets (will only be correct for one level of nesting - needs to be modified if we're having multiple levels in the future)
	my $bullet_open = "";
	my $bullet_close = "";
	if ($indent > 0) {
		$bullet_open = "<ul style=\"margin:0px\"><li style=\"margin:0px\">";
		$bullet_close = "</li></ul>";
	}
	
	# Print the set attributes
	print "\t<tr$rowclass>\n";
	print "\t\t<td>$bullet_open" . $set->name() . "$bullet_close</td>\n";
	print "\t\t<td>" . $set->short_name() . "</td>\n";
	print "\t\t<td>" . $set->description() . "</td>\n";
	print "\t</tr>\n";
	
	# Get the subsets that have the current set as immediate parent
	my $subsets = $set->get_all_sub_VariationSets(1);
	
	# Call the print subroutine for each of the subsets with an increased indentation
	my $ssets;
	foreach my $sub_vs (@{$subsets}) {
		$ssets->{$sub_vs->name} = $sub_vs;
	}
	foreach my $sset_name (sort(keys(%$ssets))) {
		print_set($ssets->{$sset_name},$rowcount,$indent+1);
	}
}

sub usage {
    
    print STDOUT qq{
Usage:

  $0 SPECIES DB_HOST DB_VERSION
  
Description:

  Prints html code for a table containing the available variation sets for a species. The species 
  has to be specified on the command line as the first argument and the host database has to be 
  specified as the second argument
         
};
    
    exit(0);
}
