#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2019] EMBL-European Bioinformatics Institute
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


# Script to update the documentation page "Data description".

use strict;
use warnings;
use DBI;
use Getopt::Long;

# Print the usage instructions if run without parameters
usage() unless (scalar(@ARGV));

my ($version,$input_dir,$output_dir,$help,$host,$hlist,$user,$port,$phost,$ohost,$species);

GetOptions(
  'v=i'         => \$version,
  'i=s'         => \$input_dir,
  'o=s'         => \$output_dir,
  'help!'       => \$help,
  'host=s'      => \$host,
  'hlist=s'     => \$hlist,
  'user=s'      => \$user,
  'port=i'      => \$port,
  'phost=s'     => \$phost,
  'ohost=s'     => \$ohost,
  'species|s=s' => \$species
);

usage("input and output directories must be specified") unless ($input_dir && $output_dir);
usage("Host, port and version must be specified") unless ($host && $port && $version);
usage("Hosts list, user must be specified") unless ($hlist && $user);
usage("Previous host must be specified") unless ($phost);
usage("Host providing the ontology database must be specified") unless ($ohost);

$species ||= 'Homo_sapiens';

my $section;
my ($tmp_file, $tmp_section, $file_name);
my ($content_before, $new_content, $content_after);

my @ontologies = ('efo', 'hpo');


#### Generates the "List of species" table documentation

# Settings
$section = 'sources';
$tmp_file    = "data_desc_$section.html";
$tmp_section = "$section\_tmp.html";
$file_name   = "species_data_types.html";

print STDOUT "# Start species list...\n";
`cp $input_dir/species/$file_name $tmp_file`;
$content_before = get_content($section,'start');
$content_after  = get_content($section,'end');
`perl species_list.pl -v $version -o $tmp_section -hlist $hlist -user $user -phost $phost`;
$new_content = `cat $tmp_section`;
`rm -f $tmp_section`;
print_into_tmp_file($tmp_file,$content_before,$new_content,$content_after);

`cp $tmp_file $output_dir/$file_name`;
print STDOUT "> Species list finished\n";


#### Generates the "Variation classes" table documentation

# Settings
$section = 'classes';
$tmp_file    = "data_desc_$section.html";
$tmp_section = "$section\_tmp.html";
$file_name   = "classification.html";

print STDOUT "# Start variant classes ...\n";
`cp $input_dir/prediction/$file_name $tmp_file`;
$content_before = get_content($section,'start');
$content_after  = get_content($section,'end');
`perl generate_classes_table.pl -v $version -o $tmp_section -host $host -port $port -ohost $ohost -species $species`;
$new_content = `cat $tmp_section`;
`rm -f $tmp_section`;
print_into_tmp_file($tmp_file,$content_before,$new_content,$content_after);

`cp $tmp_file $output_dir/$file_name`;
print STDOUT "> Variant classes finished\n";


#### Generates the "Populations" table documentation

# Settings
$section = 'populations';
$tmp_file    = "data_desc_$section.html";
$tmp_section = "$section\_tmp.html";
$file_name   = "populations.html";

print STDOUT "# Start populations ...\n";
`cp $input_dir/species/$file_name $tmp_file`;
$content_before = get_content($section,'start');
$content_after  = get_content($section,'end');
`perl generate_population_table.pl -v $version -o $tmp_section -hlist $hlist -user $user`;
$new_content = `cat $tmp_section`;
`rm -f $tmp_section`;
print_into_tmp_file($tmp_file,$content_before,$new_content,$content_after);

`cp $tmp_file $output_dir/$file_name`;
print STDOUT "> Populations finished\n";


#### Generates the "Variation sets" table documentation

# Settings
$section = 'variation_sets';
$tmp_file    = "data_desc_$section.html";
$tmp_section = "$section\_tmp.html";
$file_name   = "sets.html";

print STDOUT "# Start variant sets ...\n";
`cp $input_dir/species/$file_name $tmp_file`;
$content_before = get_content($section,'start');
$content_after  = get_content($section,'end');
`perl generate_variation_set_table.pl -v $version -o $tmp_section -hlist $hlist -user $user`;
$new_content = `cat $tmp_section`;
`rm -f $tmp_section`;
print_into_tmp_file($tmp_file,$content_before,$new_content,$content_after);

`cp $tmp_file $output_dir/$file_name`;
print STDOUT "> Variant sets finished\n";


#### Generates the "Clinical significance" tables documentation

# Settings
$section = 'clin_significance';
$tmp_file    = "data_desc_$section.html";
$tmp_section = "$section\_tmp.html";
$file_name   = "phenotype_annotation.html";

print STDOUT "# Start clinical significance ...\n";
`cp $input_dir/phenotype/$file_name $tmp_file`;
$content_before = get_content($section,'start');
$content_after  = get_content($section,'end');
`perl generate_clin_significance_tables.pl -v $version -o $tmp_section -host $host -port $port -species $species`;
$new_content = `cat $tmp_section`;
`rm -f $tmp_section`;
print_into_tmp_file($tmp_file,$content_before,$new_content,$content_after);

`cp $tmp_file $output_dir/$file_name`;
print STDOUT "> Clinical significance finished\n";


#### Update the ontology versions

# Settings
$tmp_file    = "data_desc_ontology.html";
$file_name   = "phenotype_annotation.html";

print STDOUT "# Start phenotype ontology ...\n";
# Use output file from the clinical significance update
`cp $output_dir/$file_name $tmp_file`;

my $sql_onto = qq{SELECT data_version FROM ontology WHERE name=? LIMIT 1};
my $tmp_file_content = `cat $tmp_file`;
foreach my $onto (@ontologies) {
  my $sth = get_connection_and_query("ensembl_ontology_$version", $ohost, $sql_onto, [$onto]);
  my $o_version = ($sth->fetchrow_array)[0];
  $o_version =~ s/releases\///;
  if ($o_version =~ /^(\d+):(\d+):(\d+)\s/) {
    $o_version = "$3-$2-$1";
  }
  
  $tmp_file_content =~ s/<span id="$onto\_version">(\d+(-|\.)?)+<\/span>/<span id="$onto\_version">$o_version<\/span>/i;  
}
print_into_tmp_file($tmp_file,$tmp_file_content,'','');

`cp $tmp_file $output_dir/$file_name`;
print STDOUT "> Phenotype ontology finished\n";

print STDOUT "\n>>> End of script\n";


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
    -i              Path to the input directory (Required)
    -o              Path to the output directory (Required)
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
