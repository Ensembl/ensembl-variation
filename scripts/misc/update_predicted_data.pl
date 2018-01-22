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
use Bio::EnsEMBL::Registry;
use Getopt::Long;

# Print the usage instructions if run without parameters
usage() unless (scalar(@ARGV));

my ($version,$input_file,$output_file,$help,$hlist,$user,$port,$pswd,$species,$web_colour_file,$web_mapping_colour);

GetOptions(
  'v=i'            => \$version,
  'i=s'            => \$input_file,
  'o=s'            => \$output_file,
  'help!'          => \$help,
  'hlist=s'        => \$hlist,
  'user=s'         => \$user,
  'port=i'         => \$port,
  'species|s=s'    => \$species,
#  'colour_file=s'  => \$web_colour_file,
#  'mapping_file=s' => \$web_mapping_colour
);

usage("input and output files must be specified") unless ($input_file && $output_file);
usage("Hosts list, user must be specified") unless ($hlist && $user && $version);

$species ||= 'Homo_sapiens';
$port    ||= 3306;
$pswd    ||= '';
my @hostnames = split /,/, $hlist;
my $tmp_file    = 'predicted_data_tmp.html';
my $tmp_section = 'section_tmp.html';
`cp $input_file $tmp_file`;


my $section;
my ($content_before, $new_content, $content_after);

# Generates the "List of consequences" table documentation
#$section = 'consequences';
#$content_before = get_content($section,'start');
#$content_after  = get_content($section,'end');
#`perl generate_consequence_table.pl -o $tmp_section -colour_file $web_colour_file -mapping_file $web_mapping_colour`;
#$new_content = `cat $tmp_section`;
#`rm -f $tmp_section`;
#print_into_tmp_file($tmp_file,$content_before,$new_content,$content_after);


## SIFT and PolyPhen ##

my @versions = ('sift_version', 'sift_protein_db_version', 'polyphen_version', 'polyphen_release');
my %tool_versions;
my %sift_species;
my %polyphen_species;


my $sql  = qq{SHOW DATABASES LIKE '\%variation\_$version%'};
my $sql2 = qq{SELECT meta_key,meta_value FROM meta WHERE meta_key IN ('}.join("','",@versions).qq{')};
my $sql3 = qq{SELECT meta_value FROM meta WHERE meta_key=?};

foreach my $hostname (@hostnames) {
  my $database = "";
  my $sth = get_connection_and_query($database, $hostname, $sql);
  
  # loop over databases
  while (my ($dbname) = $sth->fetchrow_array) {
    next if ($dbname =~ /^master_schema/);
    next if ($dbname =~ /sample$/);
    
    print "$dbname\n";
    $database = $dbname;
    $dbname =~ /^(.+)_variation/;
    my $s_name = $1;
    
    if ($s_name eq lc($species)) {
      my $sth2 = get_connection_and_query($database, $hostname, $sql2);
      while (my ($key,$value) = $sth2->fetchrow_array) {
        $tool_versions{$key} = $value;
      }
      $sth2->finish();
    }
    
    # SIFT
    my $sth_sift = get_connection_and_query($database, $hostname, $sql3, [$versions[0]]);
    if ($sth_sift->fetchrow_array) {
      $sift_species{$s_name} = 1;
    }
    $sth_sift->finish();
    
    # PolyPhen
    my $sth_polyphen = get_connection_and_query($database, $hostname, $sql3, [$versions[2]]);
    if ($sth_polyphen->fetchrow_array) {
      $polyphen_species{$s_name} = 1;
    }
    $sth_polyphen->finish();
  }
}  
    

# Update SIFT
$section = 'sift_version';
if ($tool_versions{$section}) {
  $content_before = get_content($section,'start');
  $content_after  = get_content($section,'end');
  my $sift_version = $tool_versions{$section};
     $sift_version =~ s/sift//;
  $new_content = $sift_version;
  print_into_tmp_file($tmp_file,$content_before,$new_content,$content_after);
}

$section = 'sift_protein_db_version';
if ($tool_versions{$section}) {
  $content_before = get_content($section,'start');
  $content_after  = get_content($section,'end');
  my $sift_pr_version = $tool_versions{$section};
  if ($sift_pr_version =~ /UniRef90/) {
    $sift_pr_version =~ s/UniRef90/UniRef90 (release/;
    $sift_pr_version .= ')';
  }    
  $new_content = $sift_pr_version;
  print_into_tmp_file($tmp_file,$content_before,$new_content,$content_after);
}

# Update SIFT species list
$section = 'sift_species';
if (scalar(%sift_species)) {
  $content_before = get_content($section,'start');
  $content_after  = get_content($section,'end');
  $new_content = print_list_of_species(\%sift_species);
  print_into_tmp_file($tmp_file,$content_before,$new_content,$content_after);
}

# Update PolyPhen
$section = 'polyphen_version';
if ($tool_versions{$section}) {
  $content_before = get_content($section,'start');
  $content_after  = get_content($section,'end');
  my $polyphen_version = $tool_versions{$section};
  $new_content = $polyphen_version;
  print_into_tmp_file($tmp_file,$content_before,$new_content,$content_after);
}

# Update PolyPhen
$section = 'polyphen_release';
if ($tool_versions{$section}) {
  $content_before = get_content($section,'start');
  $content_after  = get_content($section,'end');
  my $polyphen_release = $tool_versions{$section};
  $new_content = $polyphen_release;
  print_into_tmp_file($tmp_file,$content_before,$new_content,$content_after);
}

# Update PolyPhen species list
$section = 'polyphen_species';
if (scalar(%polyphen_species)) {
  $content_before = get_content($section,'start');
  $content_after  = get_content($section,'end');
  $new_content = print_list_of_species(\%polyphen_species);
  print_into_tmp_file($tmp_file,$content_before,$new_content,$content_after);
}

`cp $tmp_file $output_file`;
`rm -f $tmp_file`;


sub get_content {
  my $section = shift;
  my $type    = shift;
  
  my $anchor = "<!-- $section - $type -->";
  
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

sub print_list_of_species {
  my $spe_list = shift;
  my $max_row = 5;
  
  my $header = qq{
<div style="float:left;font-style:italic">
  <ul style="margin-bottom:0px">};
  my $html = $header;
  
  my $count_row = 1;
  foreach my $species (sort(keys(%$spe_list))) {
    $species =~ s/_/ /g;
    $species = ucfirst($species);
    if ($count_row == $max_row) {
      $html .= qq{  </ul>\n</div>$header};
      $count_row = 1;
    }
    $html .= ($count_row == 1) ? qq{    <li style="margin-top:0px">} : '    <li>';
    $html .= qq{$species</li>};
    $count_row ++;
  }
  $html .= qq{  </ul>\n</div>};
  return $html;
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
  my $dsn = "DBI:mysql:$dbname:$host:$port";
  my $dbh = DBI->connect($dsn, $user, $pswd) or die "Connection failed";

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
    -species        Species name. 'Homo_sapiens' by default (optional)
    -hlist          The list of host names where the new databases are stored, separated by a coma,
                    e.g. ensembldb.ensembl.org1, ensembldb.ensembl.org2 (Required)
    -user           MySQL user name (Required)
    -port           MySQL port. 3306 by default (optional)
  } . "\n";
#    -colour_file    If you want to use directly the colours from the web colours configuration file
#                    instead of the almost-up-to-date-colour-hash \%colour hash. (optional)
#                    Usually, you can find the colour configuration file in:
#                    ensembl-webcode/conf/ini-files/COLOUR.ini 
#    -mapping_file   Web module to map the colour names to the corresponding hexadecimal code. (optional)
#                    Useful because some colour names are internal to Ensembl and won't be displayed in the 
#                    documentation pages (i.e. not using the perl modules).
#                    The module ColourMap.pm can be find in:
#                    ensembl-webcode/modules/Sanger/Graphics/ColourMap.pm
#  } . "\n";
  exit(0);
}
