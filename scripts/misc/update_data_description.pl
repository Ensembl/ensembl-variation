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


# Script to update the documentation page "Data description".

use strict;
use warnings;
use DBI;
use Getopt::Long;

# Print the usage instructions if run without parameters
usage() unless (scalar(@ARGV));

my ($version,$input_file,$output_file,$help,$host,$hlist,$user,$port,$phost,$ohost,$species);

GetOptions(
  'v=i'         => \$version,
  'i=s'         => \$input_file,
  'o=s'         => \$output_file,
  'help!'       => \$help,
  'host=s'      => \$host,
  'hlist=s'     => \$hlist,
  'user=s'      => \$user,
  'port=i'      => \$port,
  'phost=s'     => \$phost,
  'ohost=s'     => \$ohost,
  'species|s=s' => \$species
);

usage("input and output files must be specified") unless ($input_file && $output_file);
usage("Host, port and version must be specified") unless ($host && $port && $version);
usage("Hosts list, user must be specified") unless ($hlist && $user);
usage("Previous host must be specified") unless ($phost);
usage("Host providing the ontology database must be specified") unless ($ohost);

$species ||= 'Homo_sapiens';
my $tmp_file    = 'data_desc_tmp.html';
my $tmp_section = 'section_tmp.html';
`cp $input_file $tmp_file`;

my @ontologies = ('efo', 'hpo');

my $section;
my ($content_before, $new_content, $content_after);

# Generates the "List of species" table documentation
$section = 'sources';
$content_before = get_content($section,'start');
$content_after  = get_content($section,'end');
`perl species_list.pl -v $version -o $tmp_section -hlist $hlist -user $user -phost $phost`;
$new_content = `cat $tmp_section`;
`rm -f $tmp_section`;
print_into_tmp_file($tmp_file,$content_before,$new_content,$content_after);


# Generates the "Variation classes" table documentation
$section = 'classes';
$content_before = get_content($section,'start');
$content_after  = get_content($section,'end');
`perl generate_classes_table.pl -v $version -o $tmp_section -host $host -port $port -ohost $ohost -species $species`;
$new_content = `cat $tmp_section`;
`rm -f $tmp_section`;
print_into_tmp_file($tmp_file,$content_before,$new_content,$content_after);


# Generates the "Populations" table documentation
$section = 'populations';
$content_before = get_content($section,'start');
$content_after  = get_content($section,'end');
`perl generate_population_table.pl -v $version -o $tmp_section -hlist $hlist -user $user`;
$new_content = `cat $tmp_section`;
`rm -f $tmp_section`;
print_into_tmp_file($tmp_file,$content_before,$new_content,$content_after);


# Generates the "Variation sets" table documentation
$section = 'variation_sets';
$content_before = get_content($section,'start');
$content_after  = get_content($section,'end');
`perl generate_variation_set_table.pl -v $version -o $tmp_section -host $host -port $port -species $species`;
$new_content = `cat $tmp_section`;
`rm -f $tmp_section`;
print_into_tmp_file($tmp_file,$content_before,$new_content,$content_after);


# Generates the "Clinical significance" tables documentation
$section = 'clin_significance';
$content_before = get_content($section,'start');
$content_after  = get_content($section,'end');
`perl generate_clin_significance_tables.pl -v $version -o $tmp_section -host $host -port $port -species $species`;
$new_content = `cat $tmp_section`;
`rm -f $tmp_section`;
print_into_tmp_file($tmp_file,$content_before,$new_content,$content_after);


# Update the ontology versions
my $sql_onto = qq{SELECT data_version FROM ontology WHERE name=? LIMIT 1};
my $tmp_file_content = `cat $tmp_file`;
foreach my $onto (@ontologies) {
  my $sth = get_connection_and_query("ensembl_ontology_$version", $ohost, $sql_onto, [$onto]);
  my $o_version = ($sth->fetchrow_array)[0];
  $o_version =~ s/releases\///;
  
  $tmp_file_content =~ s/<span id="$onto\_version">(\d+(-|\.)?)+<\/span>/<span id="$onto\_version">$o_version<\/span>/i;  
}
print_into_tmp_file($tmp_file,$tmp_file_content,'','');


`cp $tmp_file $output_file`;



#---------#
# Methods #
#---------#

sub get_content {
  my $section = shift;
  my $type    = shift;
  
  my $anchor = "<!-- Data $section - $type -->";
  
  my $line = `grep -m1 -n '$anchor' $tmp_file`;
  die "Can't find the anchor '$anchor' in the file" if (!$line || $line eq '');
  $line =~ /^(\d+):/;
  my $line_number = $1;
  my $content;  
  if ($type eq 'start') {
    $content = `head -n$line_number $tmp_file`;
    if ($content !~ /$anchor(\n?)$/) {
      $content = (split("$anchor", $content))[0].$anchor;
    }
  }
  else {
    my $lines_count = (split(' ',`wc -l $tmp_file`))[0];
    $line_number = $lines_count - $line_number + 1;
    $content = `tail -n$line_number $tmp_file`;
    if ($content !~ /^$anchor/) {
      $content = $anchor.(split("$anchor", $content))[1];
    } 
  }
  return $content;
}

sub print_into_tmp_file {
  my $tmp    = shift;
  my $before = shift;
  my $new    = shift;
  my $after  = shift;

  open  TMP, "> $tmp" or die $!;
  print TMP  $before;
  print TMP  $new;
  print TMP  $after;
}

# Connects and execute a query
sub get_connection_and_query {
  my $dbname = shift;
  my $host  = shift;
  my $sql    = shift;
  my $params = shift;

  # DBI connection 
  my $dsn = "DBI:mysql:$dbname:$host";
  my $dbh = DBI->connect($dsn, $user, '') or die "Connection failed";

  my $sth = $dbh->prepare($sql);
  if ($params) {
    $sth->execute(join(',',@$params));
  }
  else {
    $sth->execute;
  }
  return $sth;
}

sub usage {
  my $msg = shift;
  print qq{
  $msg
  Usage: perl update_data_description.pl [OPTION]
  
  Update the page "data_description.html" (under public-plugins/ensembl/htdocs/info/genome/variation/).
  
  Options:

    -help           Print this message
      
    -v              Ensembl version, e.g. 65 (Required)
    -i              Path to the data_description.html file (Required)
    -o              An HTML output file name (Required)
    -host           Host of the human database (Required)
    -port           MySQL port of the human database (Required)
    -species        Species name. 'Homo_sapiens' by default (optional)
    -hlist          The list of host names where the new databases are stored, separated by a coma,
                    e.g. ensembldb.ensembl.org1:1234, ensembldb.ensembl.org2:1234 (Required)
    -phost          Host name where the previous databases are stored, e.g. ensembldb.ensembl.org  (Required)
    -ohost          Host name where the ontology database is stored, with the port, e.g. ensembldb.ensembl.org:1234 (Required)
    -user           MySQL user name (Required)
  } . "\n";
  exit(0);
}
