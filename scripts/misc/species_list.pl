# Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

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
my ($e_version,$html_file,$hlist,$user,$port,$help);
## EG options
my ($site, $etype);

usage() if (!scalar(@ARGV));
 
GetOptions(
     'v=s'     => \$e_version,
     'o=s'     => \$html_file,
     'help!'   => \$help,
     'hlist=s' => \$hlist,
     'user=s'  => \$user,
     'port=i'  => \$port,
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
if (!$hlist) {
  print "> Error! Please give the list of host names where the new databases are stored using the option '-hlist'\n";
  usage();
}
if (!$user) {
  print "> Error! Please give user name using the option '-user'\n";
  usage();
}

usage() if ($help);

my $server_name = 'http://static.ensembl.org';
my $ecaption = 'Ensembl';
my @hostnames = split /,/, $hlist;

if ($site) {
  $server_name = $site;
}

# Settings
my $database = "";
my $pswd = "";
my $db_type = 'variation';
my $default_port = 3306;
$port ||= $default_port;

my $html;


my %tables = ( 'Genotype - Individual' => { 'order' => 2 , 'table' => 'compressed_genotype_var'},
               'Genotype - Population' => { 'order' => 3 , 'table' => 'population_genotype'},
               'Phenotype'             => { 'order' => 4 , 'table' => 'phenotype_feature'},
               'Citation'              => { 'order' => 5 , 'table' => 'variation_citation'},
               'Structural variant'    => { 'order' => 1 , 'table' => 'structural_variation'}
             );
my %columns = ( 'SIFT'     => {'order' => 1 ,'table' => 'meta', 'column' => 'meta_key', 'value' => 'sift_version'},
                'PolyPhen' => {'order' => 2 ,'table' => 'meta', 'column' => 'meta_key', 'value' => 'polyphen_version'}
              );            
my %species_list;

my $sql  = qq{SHOW DATABASES LIKE '%$db_type\_$e_version%'};
my $sql2 = qq{SELECT count(variation_id) FROM variation};

foreach my $hostname (@hostnames) {
  
  my $sth = get_connection_and_query($database, $hostname, $sql);

  # loop over databases
  while (my ($dbname) = $sth->fetchrow_array) {
    next if ($dbname =~ /^master_schema/);
    next if ($dbname =~ /sample$/);
    
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
    
    # Count the number of variations
    my $sth2 = get_connection_and_query($dbname, $hostname, $sql2);
    my $count_var = $sth2->fetchrow_array;
    $sth2->finish;
    $species_list{$s_name}{'count'} = round_count($count_var);
  }
}

my $count_species = scalar(keys(%species_list));

# Get the populated tables by species
my $species_data_tables = get_species_data_tables();

my $species_data_columns = get_species_data_columns();

my $data_tables_header  = join("</th><th>", (sort { $tables{$a}{'order'} <=> $tables{$b}{'order'} } keys(%tables)));
my $data_columns_header = join("</th><th>", (sort { $columns{$a}{'order'} <=> $columns{$b}{'order'} } keys(%columns)));

my $html_content = qq{<table class="ss" style="width:80%"><tr class="ss_header"><th>Species</th><th>Sequence variant count</th>
                      <th>$data_tables_header</th><th>$data_columns_header</th></tr>
                     };
my $bg = '';

foreach my $sp (sort keys(%species_list)) {

  my $label = $species_list{$sp}{label};
  my $uc_sp = ucfirst($sp);      
  my $img_src = "/i/species/48/$uc_sp.png";
  my $var_count = $species_list{$sp}{'count'};
  
  $html_content .= qq{
  <tr$bg>
    <td>
      <div>
        <div style="float:left;vertical-align:middle;margin-right:4px">
          <a rel="external" href="/$uc_sp/Info/Index" title="$label Ensembl Home page" style="vertical-align:middle">
            <img src="$img_src" alt="$label" class="sp-thumb" style="vertical-align:middle;width:32px;height:32px" />
          </a>
        </div>
        <div style="float:left;margin-top:0px">
          <div class="bigtext" style="font-style:italic;margin-bottom:2px">$label</div>
          <div><a href="sources_documentation.html#$sp" style="text-decoration:none" title="$label sources list">[sources]</a></div>
        </div>
        <div style="clear:both"></div>
      </div>
    </td>
    <td>$var_count variants</td>\n};
  
  # Tables
  foreach my $type (sort { $tables{$a}{'order'} <=> $tables{$b}{'order'} } keys(%tables)) {
    my $has_data = ($species_data_tables->{$sp}{$type}) ? qq{<img src="/i/16/check.png" title="Data available" />} : '-';
    $html_content .= qq{    <td style="text-align:center">$has_data</td>\n};
  }
  # SIFT, PolyPhen
  foreach my $type (sort { $columns{$a}{'order'} <=> $columns{$b}{'order'} } keys(%columns)) {
    my $has_data = ($species_data_columns->{$sp}{$type}) ?  qq{<img src="/i/16/check.png" title="Data available" />} : '-';
    $html_content .= qq{    <td style="text-align:center">$has_data</td>\n};
  }
  
  $html_content .= qq{  </tr>};
  $bg = set_bg();
}
$html_content .= qq{</table>\n};


## HTML/output file ##
open  HTML, "> $html_file" or die "Can't open $html_file : $!";
print HTML qq{<p style="padding-top:0px;margin-top:0px">There are currently <span style="font-weight:bold;font-size:1.1em;color:#000">$count_species</span> variation databases in Ensembl:</p>\n};
print HTML $html_content;
print HTML qq{<p style="padding-top:5px">The <b>full list of species</b> and their assembly versions in Ensembl is available <a href="/info/about/species.html">here</a>.</p>\n};
close(HTML);


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


sub round_count {
  my $number = shift;
  my $new_number;
  my $count;
  if ($number =~ /^(\d+)\d{6}$/) {
    $new_number = $1;
    $count  = "> $1 million";
    $count .= 's' if ($new_number != 1);
  }
  elsif ($number =~ /^(\d+)\d{3}$/) {
    $new_number = $1;
    $count = "> $1,000";
  }
  return $count;
}


# Get the list of species where the given tables are populated
sub get_species_data_tables {

  my %species_list;
  foreach my $type (keys(%tables)) {
    my $table = $tables{$type}{'table'};
    my $sql = qq{SELECT table_schema FROM information_schema.tables WHERE table_rows>=1 AND 
                 TABLE_SCHEMA like '%$db_type\_$e_version%' AND TABLE_NAME='$table'};

  
    foreach my $hostname (@hostnames) {
      my $sth = get_connection_and_query("", $hostname, $sql);

      # loop over databases
      while (my ($dbname) = $sth->fetchrow_array) {
        next if ($dbname =~ /^master_schema/);

        $dbname =~ /^(.+)_$db_type/;
        my $s_name = $1;

        $species_list{$s_name}{$type} = 1;
      }
      $sth->finish();
    }
  }
  return \%species_list;
}

# Get the list of species where the given columns are populated
sub get_species_data_columns {

  my %species_list;
  foreach my $type (keys(%columns)) {
    my $t_name = $columns{$type}{'table'};
    my $c_name = $columns{$type}{'column'};
    my $v_name = $columns{$type}{'value'};
    my $sql_col = qq{SELECT count(*) FROM $t_name WHERE $c_name="$v_name"};

    foreach my $hostname (@hostnames) {
  
      my $sth = get_connection_and_query($database, $hostname, $sql);

      # loop over databases
      while (my ($dbname) = $sth->fetchrow_array) {
        next if ($dbname =~ /^master_schema/);
        next if ($dbname =~ /sample$/);
        
        $dbname =~ /^(.+)_variation/;
        my $s_name = $1;

        my $sth_col = get_connection_and_query($dbname, $hostname, $sql_col);
        my $col_count = ($sth_col->fetchrow_array)[0];
        $sth_col->finish();
        
        $species_list{$s_name}{$type} = 1 if ($col_count != 0);
      }
      $sth->finish();
    }
  }
  return \%species_list;
}

sub set_bg {
  return ($bg eq '') ? ' class="bg2"' : '';
}

sub usage {
  
  print qq{
  Usage: perl sources2html.pl [OPTION]
  
  Put all variation sources, for each species, into an HTML document.
  
  Options:

    -help           Print this message
      
    -v              Ensembl version, e.g. 65 (Required)
    -o              An HTML output file name (Required)      
    -hlist          The list of host names where the new databases are stored, separated by a coma,
                    e.g. ensembldb.ensembl.org1, ensembldb.ensembl.org2 (Required)
    -user           MySQL user name (Required)
    -pass           MySQL password. 3306 by default (optional)
    -site           The URL of the website (optional)
    -etype          The type of Ensembl, e.g. Plant (optional)
  } . "\n";
  exit(0);
}
