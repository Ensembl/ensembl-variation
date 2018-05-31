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
  <http://www.ensembl.org/Help/Contact>.

=cut


#ÊScript to dump out a table of variation sets that can be used in the documentation

use strict;
use warnings;
use Bio::EnsEMBL::Registry;
use DBI;
use Getopt::Long;

my $registry = 'Bio::EnsEMBL::Registry';

# Print the usage instructions if run without parameters
usage() unless (scalar(@ARGV));

my ($hlist, $db_version, $output_file, $user, $help);

GetOptions(
  'v=i'     => \$db_version,
  'o=s'     => \$output_file,
  'hlist=s' => \$hlist,
  'user=s'  => \$user,
  'help!'   => \$help
);

if (!$db_version) {
  usage("> Error! Please give an Ensembl version, using the option '-v' \n");
}
if (!$output_file) {
  usage("> Error! Please give an output file using the option '-o'\n");
}
if (!$hlist) {
  usage("> Error! Please give the list of host names where the new databases are stored using the option '-hlist'\n");
}
if (!$user) {
  usage("> Error! Please give user name using the option '-user'\n");
}


my $table_header = qq{
  <tr>
    <th>Name</th>
    <th>Short name</th>
    <th>Description</th>
  </tr>
};


# Settings
my @filters = ('fail_');
my @hostnames = split /,/, $hlist;
my $database = "";
my $pswd = "";
my $db_type = 'variation';
my $img_class = "badge-48";
my $top_species = 'human';
my %display_list;
my %species_list;
my $com_sets;
my $sets;

my $sql      = qq{SHOW DATABASES LIKE '%$db_type\_$db_version%'};
my $sql_core = qq{SELECT meta_value FROM meta WHERE meta_key="species.display_name" LIMIT 1};

foreach my $hostname (@hostnames) {
  
  # Load the registry from db
  my ($host,$port) = split(':',$hostname);
  $registry->load_registry_from_db(
        -host => $host,
        -port => $port,
        -user => $user,
        -db_version => $db_version
  );
  
  my $sth = get_connection_and_query($database, $hostname, $sql);
  
  # loop over databases
  while (my ($dbname) = $sth->fetchrow_array) {
    next if ($dbname !~ /^[a-z]+_[a-z]+_variation_\d+_\d+$/i);
    next if ($dbname =~ /^master_schema/ || $dbname =~ /^homo_sapiens_variation_\d+_37$/ || $dbname =~ /private/);
    
    print $dbname;
    $dbname =~ /^(.+)_variation/;
    my $s_name = $1;
    
    my $label_name = ucfirst($s_name);
       $label_name =~ s/_/ /g;
    $species_list{$s_name}{'label'} = $label_name;
    
    # Get species display name
    my $core_dbname = $dbname;
       $core_dbname =~ s/variation/core/i;
    my $sth_core = get_connection_and_query($core_dbname, $hostname, $sql_core);
    my $display_name = $sth_core->fetchrow_array;  
       $display_name =~ s/saccharomyces/S\./i;

    # Get a VariationSetAdaptor on the human variation database
    my $vs_adaptor = $registry->get_adaptor($s_name,'variation','variationset');

    # Get all top-level variation sets
    my $top_vss = $vs_adaptor->fetch_all_top_VariationSets();

    # Loop over the top level variation sets and recursively print the subsets
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
      # Species specific set
      if (!$is_com) {
        $sets->{$s_name}{$top_vs->short_name} = $top_vs;
        $display_list{$display_name} = $s_name;
      }
    }
    print " ... done\n";
  }
}


## Print the common table headers
my $html;
$html .= qq{
  <h2 id="commom_set">Variant sets common to all species</h2>
  <div style="margin:6px 0px 30px">
    <table id=\"variation_set_table\" class=\"ss\">
    $table_header
};

foreach my $com_set_name (sort {lc $sets->{$a}->name cmp lc $sets->{$b}->name} keys(%$com_sets)) {
  my $com_rowcount = 0;
  $html .= print_set($com_sets->{$com_set_name},\$com_rowcount);
}
$html .= qq{    </table>\n  </div>\n};


foreach my $display_name (sort { $a !~ /$top_species/i cmp $b !~ /$top_species/i || $a cmp $b } keys(%display_list)) {
  my $species = $display_list{$display_name};
  ## Print the species table headers;
  my $species_label = $species_list{$species}{'label'};
  my $id_species = ucfirst($species);
  $html .= qq{
    <div style="padding-left:0px;padding-bottom:3px">
      <a href="/$id_species/Info/Index" title="$display_name Ensembl Home page" style="vertical-align:middle" target="_blank"><img src="/i/species/$id_species.png" alt="$display_name" class="$img_class" style="float:none;margin-right:4px;vertical-align:middle" /></a>
      <h2 id="$id_species" style="display:inline;color:#333">$display_name<span class="small vdoc_species_sci_name"> ($species_label)</span</h2>
    </div>
    <div style="margin:6px 0px 30px">
      <table class="ss">
      $table_header};

  my $rowcount  = 0;
  foreach my $set_name (sort {lc $sets->{$species}->{$a}->name cmp lc $sets->{$species}->{$b}->name } keys(%{$sets->{$species}})) {
    $html .= print_set($sets->{$species}->{$set_name},\$rowcount);
  }
  $html .= qq{      </table>\n    </div>\n};
  if ($display_name =~ /$top_species/i) {
    $html .= qq{
    <div style="background-color:#F0F0F0;margin:75px 0px 35px;padding:5px;border-top:2px solid #336;border-bottom:1px solid #336">
      <h2 style="display:inline;color:#000">Variant sets for the non-$top_species species</h2>
    </div>};
  }
}

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
  my $short_name = ($set->short_name())  ? $set->short_name()  : '-';
  my $set_desc   = ($set->description()) ? $set->description() : '-';
  $html_set .= "  <tr$rowclass>\n";
  $html_set .= "    <td>$bullet_open$label$bullet_close</td>\n";
  $html_set .= "    <td>$short_name</td>\n";
  $html_set .= "    <td>$set_desc</td>\n";
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


# Connects and execute a query
sub get_connection_and_query {
  my $dbname = shift;
  my $hname  = shift;
  my $sql    = shift;
  
  my ($host, $port) = split /\:/, $hname;
  
  # DBI connection 
  my $dsn = "DBI:mysql:$dbname:$host:$port";
  my $dbh = DBI->connect($dsn, $user, $pswd) or die "Connection failed";

  my $sth = $dbh->prepare($sql);
  $sth->execute;
  return $sth;
}


sub usage {
  my $msg = shift;
  print qq{
  $msg
  Usage: perl generate_variation_set_table.pl [OPTION]
  
  Update the variation set tables in "data_description.html" (under public-plugins/ensembl/htdocs/info/genome/variation/).
  
  Options:

    -help       Print this message
      
    -v          Ensembl version, e.g. 65 (Required)
    -o          An HTML output file name (Required)
    -hlist      The list of host names (with port) where the new databases are stored, separated by a coma,
                e.g. ensembldb.ensembl.org1:1234, ensembldb.ensembl.org2:1234 (Required)
    -user       MySQL user name (Required)
  } . "\n";
  exit(0);
}
