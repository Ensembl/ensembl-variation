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



=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk.org>.

=cut


#ÊScript to dump out a table of variation sets that can be used in the documentation

use strict;
use warnings;
use Bio::EnsEMBL::Registry;
use Getopt::Long;

my $registry = 'Bio::EnsEMBL::Registry';

# Print the usage instructions if run without parameters
usage() unless (scalar(@ARGV));

my ($species, $host, $port, $db_version, $output_file, $help);

GetOptions(
  'v=i'         => \$db_version,
  'o=s'         => \$output_file,
  'host=s'      => \$host,
  'port=i'      => \$port,
  'species|s=s' => \$species,
  'help!'       => \$help
);

usage ("Species, host, port, version and output_file must be specified") unless ($species && $host && $port && $db_version && $output_file);

# Filters
my @filters = ('fail_');

# Load the registry from db
$registry->load_registry_from_db(
    -host => $host,
    -port => $port,
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
my $html;
$html .= "<h4>Variant sets common to all species</h4>\n";
$html .= "<table id=\"variation_set_table\" class=\"ss\">\n";
$html .= "$table_header\n";

foreach my $com_set_name (sort {lc $sets->{$a}->name cmp lc $sets->{$b}->name} keys(%$com_sets)) {
  $html .= print_set($com_sets->{$com_set_name},\$com_rowcount);
}
$html .= "</table>\n";


## Print the human specific table headers
$html .= "<br />\n<h4>Variant sets specific to Human</h4>\n";
$html .= "<table id=\"human_variation_set_table\" class=\"ss\">\n";
$html .= $table_header;

foreach my $set_name (sort {lc $sets->{$a}->name cmp lc $sets->{$b}->name } keys(%$sets)) {
  $html .= print_set($sets->{$set_name},\$rowcount);
}
$html .= "</table>\n";


open  OUT, "> $output_file" or die $!;
print OUT $html;
close(OUT);



# We define a function that will help us recurse over the set hierarchy and print the data   
sub print_set {
  my $set = shift;
  my $rowcount = shift;
  my $indent = shift || 0;
  
  my $html_set;
  
  # Highlight even row numbers
  ${$rowcount}++;
  my $rowclass = (${$rowcount}%2 == 0 ? " class=\"bg2\"" : "");
  
  # Put a bullet next to subsets (will only be correct for one level of nesting - needs to be modified if we're having multiple levels in the future)
  my $bullet_open = "";
  my $bullet_close = "";
  my $label = $set->name();
  if ($indent > 0) {
    $bullet_open = "<ul style=\"margin:0px\"><li style=\"margin:0px\">";
    $bullet_close = "</li></ul>";
  }
  else {
    $label = "<b>$label</b>";
  }
  
  # Print the set attributes
  $html_set .= "  <tr$rowclass>\n";
  $html_set .= "    <td>$bullet_open$label$bullet_close</td>\n";
  $html_set .= "    <td>" . $set->short_name() . "</td>\n";
  $html_set .= "    <td>" . $set->description() . "</td>\n";
  $html_set .= "  </tr>\n";
  
  # Get the subsets that have the current set as immediate parent
  my $subsets = $set->get_all_sub_VariationSets(1);
  
  # Call the print subroutine for each of the subsets with an increased indentation
  my $ssets;
  foreach my $sub_vs ( sort {$a->name cmp $b->name} @{$subsets}) {
    $ssets->{$sub_vs->name} = $sub_vs;
  }
  foreach my $sset_name (sort {$a cmp $b} keys(%$ssets)) {
    $html_set .= print_set($ssets->{$sset_name},$rowcount,$indent+1);
  }
  return $html_set;
}

sub usage {
  my $msg = shift;
  print qq{
  $msg
  Usage: perl generate_variation_set_table.pl [OPTION]
  
  Update the variation set tables in "data_description.html" (under public-plugins/ensembl/htdocs/info/genome/variation/).
  
  Options:

    -help           Print this message
      
    -v              Ensembl version, e.g. 65 (Required)
    -o              An HTML output file name (Required)
    -host           Host of the human database (Required)
    -port           Human database port (Required)
    -species        Species name (Required) 
  } . "\n";
  exit(0);
}
