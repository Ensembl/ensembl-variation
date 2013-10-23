# Copyright 2013 Ensembl
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

# Script to generate an HTML page containing the variation sources of each species


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk.org>.

=cut


use Bio::EnsEMBL::Registry;
use DBI;
use strict;
use POSIX;
use Getopt::Long;

###############
### Options ###
###############
my ($e_version,$html_file,$hlist,$phost,$help);
## EG options
my ($site, $etype);

usage() if (!scalar(@ARGV));
 
GetOptions(
     'v=s'     => \$e_version,
     'o=s'     => \$html_file,
     'help!'   => \$help,
     'hlist=s' => \$hlist,
     'phost=s' => \$phost,
     'site=s'  => \$site,
		 'etype=s' => \$etype
);

if (!$e_version) {
  print "> Error! Please give an Ensembl version, using the option '-v' \n";
  usage();
}
if (!$html_file) {
  print "> Error! Please give an output file using the option '-o'\n";
  usage();
}
if (!$phost) {
  print "> Error! Please give host name where the previous databases are stored using the option '-phost'\n";
  usage();
}
if (!$hlist) {
  print "> Error! Please give the list of host names where the new databases are stored using the option '-hlist'\n";
  usage();
}

usage() if ($help);

my $server_name = 'http://static.ensembl.org';
my $ecaption = 'Ensembl';
my $previous_host = $phost;
my @hostnames = split /,/, $hlist;

if ($site) {
  $server_name = $site;
}

# Settings
my $database = "";
my $login = "ensro";
my $pswd = "";
my $sep = "\t";
my %colors = ( 'version' => '#090', 'source'  => '#00F' );
my $nb_col = 4;

my $html = 


my $html_content = '<table><tr><td>';
my %species_list;

foreach my $hostname (@hostnames) {

  my $sql = qq{SHOW DATABASES LIKE '%variation_$e_version%'};
  my $sth = get_connection_and_query($database, $hostname, $sql);

  # loop over databases
  while (my ($dbname) = $sth->fetchrow_array) {
    next if ($dbname =~ /^master_schema/);
    print $dbname;
    $dbname =~ /^(.+)_variation/;
    my $s_name = $1;
    
    if ($etype) { # EG site - need to filter out species
      my $img_thumb = sprintf qq{eg-plugins/%s/htdocs/img/species/thumb_%s.png}, $etype, ucfirst($s_name);
      #  print "- checking for $img_thumb ... ";
      if (! -e $img_thumb) {
        print "\t... skipping \n";
        next;
      } 
    }
    print "\n";
		
		my $label_name = ucfirst($s_name);
		   $label_name =~ s/_/ /g;
		$species_list{$s_name}{label} = $label_name;   
    
    # Previous database (and sources)
    my $p_version = $e_version-1;
    my $sql3 = qq{SHOW DATABASES LIKE '%$s_name\_variation_$p_version%'};
    my $sth3 = get_connection_and_query($database, $previous_host, $sql3);
    my $p_dbname = $sth3->fetchrow_array;
    
    $species_list{$s_name}{new} = 1 if (!$p_dbname);
  }
}

my $nb_sp_by_col = ceil(scalar(keys(%species_list))/$nb_col);
my $count_rows = 0;

foreach my $sp (sort keys(%species_list)) {
  if ($count_rows == $nb_sp_by_col) {
	  $html_content .= qq{\n  </td>\n  <td style="width:8px"></td>\n  <td>\n};
		$count_rows = 0;
	}
	my $label = $species_list{$sp}{label};
	my $img_src = "/i/species/48/".ucfirst($sp).".png";
	
	if ($species_list{$sp}{new}) {
	  $html_content .= qq{
    <div style="margin-bottom:5px;padding-right:4px;overflow:hidden;background-color:#336;padding:1px;color:#FFF">
	    <div style="float:left">
	      <img src="$img_src" alt="$label" class="sp-thumb" style="float:left;margin-right:4px;vertical-align:middle" />
	    </div>
	    <div style="float:left;padding-top:2px;padding-left:2px">
        <span style="font-style:italic">$label</span><br />
        <span><a href="sources_documentation.html#$sp" style="color:#AAF">[sources]</a></span><br />
        <span style="font-weight:bold;color:#F00">New species</span>
      </div>
    </div>\n};
	}
	else {
    $html_content .= qq{
    <div style="margin-bottom:5px;padding-right:4px;overflow:hidden">
	    <img src="$img_src" alt="$label" class="sp-thumb" style="float:left;margin-right:4px;vertical-align:middle" />
	    <div style="float:left;padding-top:6px">
        <div style="font-style:italic;margin-bottom:3px">$label</div>
        <span><a href="sources_documentation.html#$sp">[sources]</a></span>
      </div>
    </div>\n};
	}
	$count_rows ++;
}
$html_content .= qq{</td></tr></table>\n};


## HTML/output file ##
open  HTML, "> $html_file" or die "Can't open $html_file : $!";
print HTML $html_content;
close(HTML);


# Connects and execute a query
sub get_connection_and_query {
  my $dbname = shift;
  my $hname  = shift;
  my $sql    = shift;
  
  my ($host, $port) = split /\:/, $hname;

  # DBI connection 
  my $dsn = "DBI:mysql:$dbname:$host:$port";
  my $dbh = DBI->connect($dsn, $login, $pswd) or die "Connection failed";

  my $sth = $dbh->prepare($sql);
  $sth->execute;
  return $sth;
}


sub usage {
  
  print qq{
  Usage: perl sources2html.pl [OPTION]
  
  Put all variation sources, for each species, into an HTML document.
  
  Options:

    -help           Print this message
      
    -v              Ensembl version, e.g. 65 (Required)
    -o              An HTML output file name (Required)      
    -phost          Host name where the previous databases are stored, e.g. ensembldb.ensembl.org  (Required)
    -hlist          The list of host names where the new databases are stored, separated by a coma,
                    e.g. ensembldb.ensembl.org1, ensembldb.ensembl.org2 (Required)
    -site           The URL of the website (optional)
    -etype          The type of Ensembl, e.g. Plant (optional)
  } . "\n";
  exit(0);
}
