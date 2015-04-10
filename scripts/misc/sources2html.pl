# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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
my ($e_version,$html_file,$source_id,$source,$s_version,$s_description,$s_url,$s_type,$s_status,$s_data_types,$s_order,$hlist,$phost,$help);
my ($set_id,$set_name,$set_description);

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
  print STDERR "> Error! Please give an Ensembl version, using the option '-v' \n";
  usage();
}
if (!$html_file) {
  print STDERR "> Error! Please give an output file using the option '-o'\n";
  usage();
}
if (!$phost) {
  print STDERR "> Error! Please give host name where the previous databases are stored using the option '-phost'\n";
  usage();
}
if (!$hlist) {
  print STDERR "> Error! Please give the list of host names where the new databases are stored using the option '-hlist'\n";
  usage();
}

usage() if ($help);

my %top_species = ('homo_sapiens' => 1);
my $server_name = 'http://static.ensembl.org';
my $ecaption = 'Ensembl';
my $previous_host = $phost;
my @hostnames = split /,/, $hlist;

if ($site) {
  $server_name = $site;
}
if ($etype) {
  $ecaption .= ' '.ucfirst($etype);
}
# Settings
my $database = "";
my $login = "ensro";
my $pswd = "";
my $sep = "\t";
my $start = 0;
my %colours = ( 'version'     => '#090',
                'source'      => '#00F',
                'lot_million' => '#800',
                'few_million' => '#007',
                'thousand'    => '#006266', 
                'hundred'     => '#070'
              );

my $phen_icon = '/i/val/var_phenotype_data_small.png';
my $internal_link = '/i/16/internal_link.png';

my %data_type_example = (
  'variation'            => {
                             'sql'       => qq{SELECT name FROM variation WHERE somatic=0 AND variation_id NOT IN (SELECT variation_id FROM failed_variation) AND source_id=? LIMIT 1},
                             'sql_som'   => qq{SELECT name FROM variation WHERE variation_id NOT IN (SELECT variation_id FROM failed_variation) AND source_id=? LIMIT 1},
                             'count_spe' => qq{SELECT source_id, COUNT(variation_id) FROM variation GROUP BY source_id},
                             'url'       => 'Variation/Explore?v='
                             },
  'variation_synonym'    => {
                             'sql'       => qq{SELECT v.name FROM variation v, variation_synonym s WHERE v.variation_id=s.variation_id AND 
                                           s.source_id= ? AND v.variation_id NOT IN (SELECT variation_id FROM failed_variation) LIMIT 1},
                             'count_spe' => qq{SELECT source_id, COUNT(variation_synonym_id) FROM variation_synonym GROUP BY source_id},
                             'url'       => 'Variation/Explore?v='
                            },
  'structural_variation' => {
                             'sql'       => qq{SELECT variation_name FROM structural_variation WHERE is_evidence=0 AND somatic=0 AND structural_variation_id NOT IN 
                                           (SELECT structural_variation_id FROM failed_structural_variation) AND source_id=? LIMIT 1},
                             'count_spe' => qq{SELECT source_id, COUNT(structural_variation_id) FROM structural_variation WHERE is_evidence=0 GROUP BY source_id},            
                             'url'       => 'StructuralVariation/Explore?sv='
                            },
  'phenotype_feature'    => {
                             'sql'       => qq{SELECT object_id, type, phenotype_id FROM phenotype_feature WHERE source_id=? AND is_significant=1 AND type!="SupportingStructuralVariation" LIMIT 1},
                             'count_spe' => qq{SELECT source_id, COUNT(phenotype_feature_id) FROM phenotype_feature GROUP BY source_id},
                             'types'     => qq{SELECT source_id, GROUP_CONCAT(DISTINCT type ORDER BY type ASC SEPARATOR ', ')
                                               FROM phenotype_feature WHERE type!="SupportingStructuralVariation" GROUP BY source_id},
                             'Variation'           => 'Variation/Phenotype?v=',
                             'StructuralVariation' => 'StructuralVariation/Phenotype?v=',
                             'Gene'                => 'Gene/Phenotype?g=',
                             'QTL'                 => 'Phenotype/Locations?ph='
                            },
  'variation_set'        => {
                             'sql'   => qq{SELECT v.name FROM variation v, variation_set_variation s WHERE v.variation_id=s.variation_id AND v.variation_id NOT IN 
                                           (SELECT variation_id FROM failed_variation) AND s.variation_set_id=? LIMIT 1},
                             'count' => qq{SELECT COUNT(variation_id) FROM variation_set_variation WHERE variation_set_id=?},
                             'url'   => 'Variation/Explore?v='
                            },                         
);

##############
### Header ###
##############
my $html_header = q{
<html>
<head>
  <title>Variation Sources</title>
  <script type="text/javascript">
    window.onload = function() {
      $('.conhelp').helptip({'track': true});
    };
  </script>
</head>

<body>

<div>
};

my $html_title = qq{
  <div style="float:left;width:75%">
    <h1 style="margin-top:15px">Ensembl Variation - Sources Documentation</h1>

    <h2>List of Variation sources for each species - $ecaption $e_version</h2>  
};


##############
### Footer ###
##############
my $html_footer = qq{
  </div>
</div>
</body>
</html>};


############
### Main ###
############

my $html_top_content = '';
my $html_content = '';
my %top_species_list;
my @species_list;
my %species_news;

my $cols_sql2  = 'source_id, name, version, description, url, type, somatic_status, data_types';
my $term_sql2  = 'dbSNP';
my $condition1 = qq{SELECT $cols_sql2, 1 AS ordering FROM source WHERE name = "$term_sql2"};
my $condition2 = qq{SELECT $cols_sql2, 2 AS ordering FROM source WHERE name like "%$term_sql2%" AND name != "$term_sql2"};
my $condition3 = qq{SELECT $cols_sql2, 3 AS ordering FROM source WHERE name not like "%$term_sql2%"};

my $sql2 = qq{$condition1 UNION $condition2 UNION $condition3 ORDER BY ordering, name};
my $sql4 = qq{SELECT name, version FROM source};

my $sql2b = qq{SELECT variation_set_id, name, description FROM variation_set WHERE 
               (name like "%illumina%" OR name like "%affymetrix%" OR description like "%illumina%" OR description like "%affymetrix%")
                AND name NOT IN (SELECT name FROM source)};
my $sql4b = $sql2b;

foreach my $hostname (@hostnames) {

  my $sql = qq{SHOW DATABASES LIKE '%variation_$e_version%'};
  my $sth = get_connection_and_query($database, $hostname, $sql);
  my $db_found = 0;
  
  # loop over databases
  while (my ($dbname) = $sth->fetchrow_array) {
    next if ($dbname =~ /^master_schema/);
    $db_found ++;
    print STDERR $dbname;
    $dbname =~ /^(.+)_variation/;
    my $s_name = $1;

    if ($etype) { # EG site - need to filter out species
      my $img_thumb = sprintf qq{eg-plugins/%s/htdocs/img/species/thumb_%s.png}, $etype, ucfirst($s_name);
      if (! -e $img_thumb) {
        print STDERR "\t... skipping \n";
        next;
      } 
    }
    print STDERR "\n";
    # Get list of sources from the new databases
    my $sth2 = get_connection_and_query($dbname, $hostname, $sql2);
    $sth2->bind_columns(\$source_id,\$source,\$s_version,\$s_description,\$s_url,\$s_type,\$s_status,\$s_data_types,\$s_order);
    
    my $sth2b = get_connection_and_query($dbname, $hostname, $sql2b);
    $sth2b->bind_columns(\$set_id,\$set_name,\$set_description);
    
    
    # Previous database (and sources)
    my $p_version = $e_version-1;
    my $sql3 = qq{SHOW DATABASES LIKE '%$s_name\_variation_$p_version%'};
    my $sth3 = get_connection_and_query($database, $previous_host, $sql3);
    my $p_dbname = $sth3->fetchrow_array;
 
    my %p_list;
    my %p_set_list;
    my $is_new_species = 0;
    if ($p_dbname) {
      # Previous sources
      my $sth4 = get_connection_and_query($p_dbname, $previous_host, $sql4);
      while (my @p = $sth4->fetchrow_array) {
        $p_list{$p[0]} = $p[1];
      }
      # Previous sets
      my $sth4b = get_connection_and_query($p_dbname, $previous_host, $sql4b);
      while (my @s = $sth4b->fetchrow_array) {
        $p_set_list{$s[1]} = $s[2];
      }
    }
    else {
      $is_new_species = 1;
    }
    
    # Display the species at the top of the list
    if ($top_species{$s_name}) {
      $html_top_content .= source_table($s_name,$sth2,$sth2b,$is_new_species,$dbname,$hostname,\%p_list,\%p_set_list);
    }
    else {
      $html_content .= source_table($s_name,$sth2,$sth2b,$is_new_species,$dbname,$hostname,\%p_list,\%p_set_list);
    }
    
    $start = 1 if ($start == 0);
  }
   die ("No variation databases found on $hostname for the version $e_version!") if ($db_found == 0);
}

my $html_menu = create_menu();

if ($html_content ne '') {
  $html_top_content .= qq{
    <div style="background-color:#F0F0F0;margin:50px 0px 25px;padding:5px;border-top:2px solid #336;border-bottom:1px solid #336">
      <h2 style="display:inline;color:#000">Others species</h2>
    </div>
  };
}


## HTML/output file ##
open  HTML, "> $html_file" or die "Can't open $html_file : $!";
print HTML $html_header."\n";
print HTML $html_menu."\n";
print HTML $html_title."\n";
print HTML $html_top_content."\n";
print HTML $html_content."\n";
print HTML $html_footer."\n";
close(HTML);


###############
### Methods ###
###############

sub source_table {
  my $name       = shift;
  my $sth        = shift;
  my $sth_set    = shift;
  my $is_new     = shift;
  my $db_name    = shift;
  my $hostname   = shift;
  my $p_list     = shift;
  my $p_set_list = shift;

  
  my $s_name = ucfirst($name);
  my $species = $s_name;
     $species =~ s/_/ /;
  my $s_name_id = $name;
  
  if ($top_species{$name}) {
    $top_species_list{$name} = {name => $species, s_name => $s_name, anchor => $s_name_id};
  }
  else {
    push (@species_list,{name => $species, s_name => $s_name, anchor => $s_name_id});
  }
  
  my $html = qq{<!-- $species -->};
  if ($is_new) {
    $html .= qq{
    <div style="padding-left:0px;padding-bottom:1px">
      <a rel="external" href="/$s_name/Info/Index" title="$species Ensembl Home page" style="vertical-align:middle"><img src="/i/species/48/$s_name.png" alt="$species" class="sp-thumb" style="float:none;margin-right:4px;padding:2px;vertical-align:middle;background-color:#00F" /></a>
      <h2 id="$s_name_id" style="display:inline;color:#333">$species</h2><span style="padding-left:20px;color:#00F;font-weight:bold">New species!</span>
    </div>
    };
  }
  else {
    $html .= qq{
    <div style="padding-left:0px;padding-bottom:3px">
      <a rel="external" href="/$s_name/Info/Index" title="$species Ensembl Home page" style="vertical-align:middle"><img src="/i/species/48/$s_name.png" alt="$species" class="sp-thumb" style="float:none;margin-right:4px;vertical-align:middle" /></a>
      <h2 id="$s_name_id" style="display:inline;color:#333">$species</h2>
    </div>
    };
  }
  
  my $source_table;
  
  my $bg = 1;
  my @p_sources = keys(%{$p_list});
  my @p_sets = keys(%{$p_set_list});
  my %other_flag;
  
  # Chip headers
  my $cbg = 1;
  my $chip_table;
     
  # LSDB headers
  my $lbg = 1;
  my $lsdb_table;
  
  
  my $counts_species = get_species_count(\%data_type_example, $s_name, $db_name, $hostname);
  my $phe_types      = get_phenotype_types(\%data_type_example, $s_name, $db_name, $hostname);

  my $type_style  = qq{style="float:left;width:65px;text-align:left"};
  my $count_style = qq{style="float:left;width:60px;text-align:right;margin:0px 5px"};
  my $eg_style    = qq{style="float:left;width:20px;text-align:center"};
  my $spaces      = "          ";
  my $phe_title   = "Provides phenotype association data";
  
  ########## Sources ##########    
  
  while ($sth->fetch) {
  
    # Check if the source or its version is new
    my $s_new      = '';
    my $s_new_type = '';
    my $s_header   = '<td style="width:4px;padding:0px;margin:0px';
    if (!grep {$_ eq $source} @p_sources) {
      $s_new_type = 'source';
    }
    elsif ($p_list->{$source} ne $s_version){
      $s_new_type = 'version';
    }
   
    if ($s_new_type) {
      $species_news{$species}{$s_new_type} = 1;
      my $borders = ";border-top:1px solid #FFF;border-bottom:1px solid #FFF";
      if ($s_type eq 'chip') {
        $s_header .= $borders if ($cbg == 1);
      }
      elsif ($s_type eq 'lsdb') {
        $s_header .= $borders if ($lbg == 1);
      }
      else {
        $s_header .= $borders if ($bg == 1);
      }
      $s_new = '<span style="color:'.$colours{$s_new_type}.'">New '.$s_new_type.'</span>' if ($s_new_type);
      $s_header .= ';background-color:'.$colours{$s_new_type};
    }
      
    $s_header .= '"></td>';

    my $sql5 = qq{SELECT variation_set_id FROM variation_set WHERE name="$source" AND name IN (SELECT name FROM source) AND 
                  (name like "%illumina%" OR name like "%affymetrix%" OR description like "%illumina%" OR description like "%affymetrix%")};

    # Source
    if ($s_url) {
      $source = qq{<a href="$s_url" style="text-decoration:none" target="_blank">$source</a>};
    }
    
    # Version
    $s_version = format_version($s_version);  
    
    # Somatic status
    my $s_somatic_status = somatic_status($s_status);
    my $is_somatic = ($s_status eq 'somatic') ? 1 : undef;
    
    # Phenotype
    my $s_phenotype = '';
    
    # Data types
    my @data_types = split(",", $s_data_types);
    my $data_type_string = '';
    my $counts;
    my $examples;
    
    foreach my $dt (@data_types) {
      $data_type_string .= qq{$spaces<div>};
      
      my $data_type_label = ucfirst($dt);
      $data_type_label =~ s/_/ /g;

      if ($dt eq 'phenotype_feature') {
        my $dt_phe_title = ($phe_types->{$dt}{$source_id}) ? "Provides ".$phe_types->{$dt}{$source_id}." phenotype association data" : $phe_title;
        $data_type_string .= qq{\n$spaces  <div $type_style><span class="_ht conhelp" title="$dt_phe_title">Phenotype</span></div>};
        $s_phenotype = qq{<img src="$phen_icon" style="border-radius:5px;border:1px solid #000" alt="$dt_phe_title" title="$dt_phe_title" />};
      }
      elsif ($dt eq 'study') {
        $data_type_string .= qq{\n$spaces  <div $type_style><span class="_ht conhelp" title="Data are grouped by study/publication">$data_type_label</span></div>};
      }
      elsif ($dt eq 'variation_synonym') {
        $data_type_string .= qq{\n$spaces  <div $type_style><span class="_ht conhelp" title="$data_type_label - Some/all variants already exist in an other source, or are redundant in this source, with different IDs">Synonym</span></div>};
      }
      elsif ($dt eq 'structural_variation') {
        $data_type_string .= qq{\n$spaces  <div $type_style><span class="_ht conhelp" title="$data_type_label">SV</span></div>};
      }
      else {
        $data_type_string .= qq{\n$spaces  <div $type_style>$data_type_label</div>};
      }
      
      
      # Count
      my $count = $counts_species->{$dt}{$source_id};
      $data_type_string .= qq{\n$spaces  <div $count_style>$count</div>};
      
      # Example
      my $somatic_example = ($is_somatic && $dt eq 'variation') ? 1 : undef;
      my $example = get_example($dt, $source_id, $s_name, $db_name, $hostname, $somatic_example);
      $data_type_string .= qq{\n$spaces  <div $eg_style>$example</div>};
      
      $data_type_string .= qq{\n$spaces  <div style="clear:both"></div>\n$spaces</div>};
    }
    
    my $sth5 = get_connection_and_query( $db_name, $hostname, $sql5);
    my ($source_var_set_id) = $sth5->fetchrow_array;
    
    # Variation source (i.e. chip data) having data in variation set as well
    if ($source_var_set_id) { # Also in variation set 
      $data_type_string .= qq{\n$spaces<div>\n$spaces  <div $type_style><span class="_ht conhelp" title="Variation set - Existing variants from 1 or several sources have been associated with this variation set">Set</span></div>};
    
      # Count
      my $count = get_species_set_count($source_var_set_id, $s_name, $db_name, $hostname);
      $data_type_string .= qq{\n$spaces  <div $count_style>$count</div>};;
   
      # Example
      my $example = get_example('variation_set', $source_var_set_id, $s_name, $db_name, $hostname);
      $data_type_string .= qq{\n$spaces  <div $eg_style>$example</div>};
    
      $data_type_string .= qq{\n$spaces  <div style="clear:both"></div>\n$spaces</div>};
    }
    
    $data_type_string = '-' if ($data_type_string eq '');
    
    $s_type = 'main' if (!defined($s_type));
    $other_flag{$s_type} = 1 if ($s_phenotype ne '' || $s_somatic_status ne '-');
    
    my $row = set_row($s_header,$source,$s_version,$s_description,$data_type_string,$s_phenotype,$s_somatic_status);
    
    # Is chip ?
    if ($s_type eq 'chip') {
      $chip_table .= qq{
      <tr class="bg$cbg">
        $row
      </tr>};
      if ($cbg == 1) { $cbg = 2; }
      else { $cbg = 1; }
    }
    # Is lsdb ?
    elsif ($s_type eq 'lsdb') {
      $lsdb_table .= qq{
      <tr class="bg$lbg">
       $row
      </tr>};
      if ($lbg == 1) { $lbg = 2; }
      else { $lbg = 1; }
    }
    else {
      $source_table .= qq{
      <tr class="bg$bg">
        $row
      </tr>};
      if ($bg == 1) { $bg = 2; }
      else { $bg = 1; }
    }
  }
  
  
  ########## Sets ##########
  
  while ($sth_set->fetch) {
   
    # Check if the set is new
    my $s_new      = '';
    my $s_new_type = '';
    my $s_header   = '<td style="width:4px;padding:0px;margin:0px';
    if (!grep {$_ eq $set_name} @p_sets) {
      $s_new_type = 'source';
    }
   
    if ($s_new_type) {
      $species_news{$species}{$s_new_type} = 1;
      my $borders = ";border-top:1px solid #FFF;border-bottom:1px solid #FFF";
      if ($s_type eq 'chip') {
        $s_header .= $borders if ($cbg == 1);
      }
      elsif ($s_type eq 'lsdb') {
        $s_header .= $borders if ($lbg == 1);
      }
      else {
        $s_header .= $borders if ($bg == 1);
      }
      $s_new = '<span style="color:'.$colours{$s_new_type}.'">New '.$s_new_type.'</span>' if ($s_new_type);
      $s_header .= ';background-color:'.$colours{$s_new_type};
    }
      
    $s_header .= '"></td>';

    
    # Display
    my $source = qq{$set_name};
    
    # Version
    my $s_version = format_version('');  
    
    # Data types
    my @data_types = split(",", $s_data_types);
    my $data_type_string = qq{\n$spaces<div>\n$spaces  <div $type_style><span class="_ht conhelp" title="Variation set - Existing variants from 1 or several sources have been associated with this variation set">Set</span></div>};
    
    # Count
    my $count = get_species_set_count($set_id, $s_name, $db_name, $hostname);
    $data_type_string .= qq{\n$spaces  <div $count_style>$count</div>};;
   
    # Example
    my $example = get_example('variation_set', $set_id, $s_name, $db_name, $hostname);
    $data_type_string .= qq{\n$spaces  <div $eg_style>$example</div>};
    
    $data_type_string .= qq{\n$spaces  <div style="clear:both"></div>\n$spaces</div>};
    
    my $row = set_row($s_header,$source,$s_version,$set_description,$data_type_string,'','');

    $chip_table .= qq{
    <tr class="bg$cbg">
      $row
    </tr>};
    if ($cbg == 1) { $cbg = 2; }
    else { $cbg = 1; }

  }
 
  # Main source header
  $source_table = table_header('Source','main',\%other_flag).$source_table;
  
  # Chip header
  $chip_table = table_header('Chip Source','chip',\%other_flag).$chip_table if ($chip_table);
     
  # LSDB header
  $lsdb_table = table_header('LSDB Source','lsdb',\%other_flag).$lsdb_table if ($lsdb_table);
  
  $html .= qq{$source_table\n    </table>\n  </div>};
  $html .= qq{$chip_table\n    </table>\n  </div>} if ($chip_table);
  $html .= qq{$lsdb_table\n    </table>\n  </div>} if ($lsdb_table);
  $html .= qq{<div style="height:20px"></div>};
  
  return $html;
}


sub format_version {
  my $version = shift;
  
  # e.g. 20110513
  if ($version =~ /(20\d{2})(\d{2})(\d{2})/) {
    $version = "$3/$2/$1";
  }
  # e.g. 201105
  elsif ($version =~ /(20\d{2})(\d{2})/) {
    $version = "$2/$1";
  }
  # e.g. 110408
  elsif ($version =~ /(\d{2})(\d{2})(\d{2})/) {
    $version = "$3/$2/20$1";
  }
  # e.g. 20121 (HGMD data version)
  elsif ($version =~ /^(20\d{2})(\d)$/) {
    $version = "$1.$2";
  }
  elsif ($version eq '') {
    $version = '-';
  }
  
  return $version;
}

# Connects and execute a query
sub get_connection_and_query {
  my $dbname = shift;
  my $hname  = shift;
  my $sql    = shift;
  my $params = shift;
  
  my ($host, $port) = split /\:/, $hname;

  # DBI connection 
  my $dsn = "DBI:mysql:$dbname:$host:$port";
  my $dbh = DBI->connect($dsn, $login, $pswd) or die "Connection failed";

  my $sth = $dbh->prepare($sql);
  if ($params) {
    $sth->execute(join(',',@$params));
  }
  else {
    $sth->execute;
  }
  return $sth;
}


sub create_menu {
  
  my $label_style = "display:inline-block;height:8px;width:8px;border-radius:4px";
  my %desc = ( 
              'version' => qq{New data version(s) for this species},
              'source'  => qq{New data source(s) for this species}
             );
  
  my $html = qq{
  <!-- Right hand side menu -->
  <div style="float:right">
  <div style="margin-left:8px;margin-top:2px;padding-bottom:2px;background-color:#F9F9F9;color:#333;border-left:1px dotted #BBB;border-right:1px dotted #BBB;border-bottom:1px dotted #BBB">
    <div style="padding:5px;font-weight:bold;color:#000;background-color:#FFF;border-top:2px solid #336;border-bottom:1px solid #336;;margin-bottom:5px">
      <img src="/i/16/info.png" style="vertical-align:top" alt="info" /> 
      Species list
    </div>
  };
  foreach my $top_sp (sort { $top_species{$a} <=> $top_species{$b} } keys(%top_species)) {
    my $species = $top_species_list{$top_sp};
    $html .= menu_list($species,$label_style,\%desc) if (defined($species));
  }
  if (scalar(keys(%top_species_list))) {
    $html .= qq{<div style="border-top:1px dotted #336;height:1px;margin:2px 0px 6px"></div>};
  }
  
  foreach my $species (@species_list) {
    $html .= menu_list($species,$label_style,\%desc);
  }
  my $v_colour  = $colours{'version'};
  my $s_colour  = $colours{'source'};
  my $lm_colour = $colours{'lot_million'};
  my $fm_colour = $colours{'few_million'};
  my $t_colour  = $colours{'thousand'};
  my $h_colour  = $colours{'hundred'};
  my $v_label   = $desc{'version'};
  my $s_label   = $desc{'source'};
  
  my $legend_div_id = 'legend';

  $html .= sprintf ( qq{
    <span style="$label_style;margin-left:5px;background-color:$v_colour"></span><small> : $v_label</small>
    <br />
    <span style="$label_style;margin-left:5px;background-color:$s_colour"></span><small> : $s_label</small>
  </div>
  <br />
  <div id="$legend_div_id" style="margin-left:8px;margin-top:2px;padding-bottom:2px;background-color:#F9F9F9;color:#333;border-left:1px dotted #BBB;border-right:1px dotted #BBB;border-bottom:1px dotted #BBB">
    <!-- Legend header -->
    <div style="padding:5px;font-weight:bold;color:#000;background-color:#FFF;border-top:2px solid #336;border-bottom:1px solid #336;;margin-bottom:2px">
      <div style="float:left">
        <img src="/i/16/info.png" style="vertical-align:top" />
        Icons legend
      </div>
      <div style="float:right">
        <a href="#top" style="text-decoration:none">[Top]</a>
      </div>
      <div style="clear:both"></div>
    </div>
    <!-- Main legend -->
    <table>
      <tr>
        <td style="padding-top:6px;text-align:center">%s   </td>
        <td style="padding-top:4px"><b>New version</b> of the data<br />source in this release<br />for the species</td>
      </tr>
      <tr>
        <td style="padding-top:6px;text-align:center">%s   </td>
        <td style="padding-top:4px"><b>New data source</b> in this<br />release for the species</td>
      </tr>
      <tr>
        <td style="padding-top:6px;text-align:center;">
          <img src="$phen_icon" style="margin-left:auto;margin-right:auto;border-radius:5px;border:1px solid #000;margin-right:1px" alt="Provides phenotype data" title="Provides phenotype data"/>
        </td>
        <td style="padding-top:4px">Source which provides<br />phenotype association data</td>
      </tr>
      <tr>
        <td style="padding-top:6px;text-align:center">%s   </td>
        <td style="padding-top:4px">The source contains only<br />germline data</td>
      </tr>
      <tr>
        <td style="padding-top:6px;text-align:center">%s    </td>
        <td style="padding-top:4px">The source contains only<br />somatic data</td>
      </tr>
      <tr>
        <td style="padding-top:6px;text-align:center">%s    </td>
        <td style="padding-top:4px">The source contains both<br />germline and somatic data</td>
      </tr>
    </table>
    
    <!-- Variant and structural variant count colour legend -->
    <div style="border-top:1px dotted #BBB;margin-top:2px;padding:4px 0px 0px">
      <span style="padding-left:4px;font-weight:bold">Data types - variants count:</span>
      <table>
        <tr>
          <td style="padding-top:4px;text-align:center">
            <span style="background-color:$lm_colour;color:#FFF;border-radius:5px;padding:1px 2px 1px 20px;cursor:default"></span>
          </td>
          <td style="padding-top:4px">greater than 10 million</td>
        </tr>
        <tr>
          <td style="padding-top:4px;text-align:center">
            <span style="background-color:$fm_colour;color:#FFF;border-radius:5px;padding:1px 2px 1px 20px;cursor:default"></span>
          </td>
          <td style="padding-top:4px">from 1 million to 9.9 million</td>
        </tr>
        <tr>
          <td style="padding-top:3px;text-align:center">
            <span style="background-color:$t_colour;color:#FFF;border-radius:5px;padding:1px 2px 1px 20px;cursor:default"></span>
          </td>
          <td style="padding-top:3px">from 1,000 to 999,999</td>
        </tr>
        <tr>
          <td style="padding-top:3px;text-align:center">
            <span style="background-color:$h_colour;color:#FFF;border-radius:5px;padding:1px 2px 1px 20px;cursor:default"></span>
          </td>
          <td style="padding-top:3px">less than 1,000</td>
        </tr>
      </table>
    </div>

    <!-- Javascript used to fix the legend on the right handside when you scroll down -->
    <script language="Javascript" type="text/javascript">
      var legend_element = document.getElementById("$legend_div_id");
      var legend_element_pos = legend_element.offsetTop;
      window.onscroll = function() {
        if (window.pageYOffset-80 > legend_element_pos)  {
          legend_element.style.position = "fixed";
          legend_element.style.top = "10px";
        }
        else {
          legend_element.style.position = "relative";
          legend_element.style.top = "auto";
        }
      };
    </script>

  </div>
</div>
  },
  qq{<div style="width:4px;height:25px;background-color:$v_colour;margin-left:9px"></div>},
  qq{<div style="width:4px;height:25px;background-color:$s_colour;margin-left:9px"></div>},
  somatic_status('germline'),
  somatic_status('somatic'),
  somatic_status('mixed'));
  return $html;
}


sub menu_list {
  my $species = shift;
  my $label_style = shift;
  my $desc = shift;
  
  my $name = $species->{name};
  my $s_name = $species->{s_name};
  my $anchor = $species->{anchor};
  my $new_data = '';
  if ($species_news{$species->{name}}) {
    my @types = sort {$b cmp $a} keys(%{$species_news{$species->{name}}});
    my $count_type = 0;
    foreach my $type (@types) {
      $count_type ++;
      my $type_margin = ($count_type == scalar(@types)) ? '2px' : '1px';
      my $label_colour = $colours{$type};
      my $label_desc = $desc->{$type};
      $new_data .= qq{<span style="$label_style;margin-left:2px;margin-right:$type_margin;background-color:$label_colour" title="$label_desc"></span>};
    }
  }
  my $img = $name;
  return qq{
  <div style="margin-left:4px;margin-bottom:5px">
    <img src="/i/species/16/$s_name.png" alt="$name" style="border-radius:4px;margin-right:4px;vertical-align:middle" />
    <a href="#$anchor" style="margin-right:2px;text-decoration:none;vertical-align:middle">$name</a>$new_data
  </div>
  };
}



sub new_source_or_version {
  my $type = shift;
  my $color = $colours{$type};
  return qq{
          <div style="color:$color;font-size:0.8em;text-align:center;margin:0px auto 0px auto;padding:0px">
            <span style="text-align:center;margin:0px;padding:0px">New</span><br /><span style="text-align:center;margin:0px;padding:0px">$type</span>
          </div>
     };
}

sub somatic_status {
  my $type = shift;
  my $html;
  if ($type eq 'germline') {
     $html .= qq{
          <div class="_ht" style="margin-left:auto;margin-right:auto;border-radius:5px;border:1px solid #000;width:20px;height:20px;background-color:#00C;" title="$type data"></div>
     };
  }
  elsif ($type eq 'somatic') {
    $html .= qq{
          <div class="_ht" style="margin-left:auto;margin-right:auto;border-radius:5px;border:1px solid #000;width:20px;height:20px;background-color:#C00;" title="$type data"></div>
    };
  }
  elsif ($type eq 'mixed') {
    $html .= qq{
          <div class="_ht" style="margin-left:auto;margin-right:auto;border-radius:5px;border:1px solid #000;width:20px;height:20px;background-color:#00C;" title="$type data">
            <div style="width:0px;height:0px;border-style:solid;border-width:0 0 20px 20px;border-color:transparent transparent #C00 transparent"></div>
          </div>
    };
  }
  else {
   return '-';
  }
}


sub table_header {
  my $name = shift;
  my $type = shift;
  my $flag = shift;
  
  my $alt_text = qq{Phenotype data, somatic/germline data, ... See the icons description on the table on the right handside of the page};
  my $header_col = qq{
    <th colspan=2 style="width:56px;text-align:center;border-left:1px solid #CCC;background-color:#BBB">
       <span class="_ht conhelp" title="$alt_text">Other</span>
    </th>};
  
  my $top_margin = ($type eq 'main') ? '6px' : '0px';
  
  my $data_type_header = qq{
     <th style="width:155px;text-align:center;border-left:1px solid #CCC;background-color:#BBB">Data type(s)
       <div>
         <div style="float:left;width:65px;text-align:center"><small>Type</small></div>
         <div style="float:left;width:70px;text-align:center"><span class="_ht conhelp" title="Variants count"><small>Count</small></span></div>
         <div style="float:left;width:20px;text-align:center"><span class="_ht conhelp" title="Example"><small>e.g.</small></span></div>
         <div style="clear:both"></div>
       </div>
     </th>
  };
  
  return qq{    <div style="margin:$top_margin 0px 30px">
      <table class="ss">
        <tr><th colspan="2">$name</th><th>Version</th><th style="max-width:800px">Description</th>$data_type_header</th>$header_col</tr>
    };
}


sub set_row {
  my $header         = shift;
  my $source         = shift;
  my $version        = shift;
  my $desc           = shift;
  my $data_type      = shift;
  my $phenotype      = shift;
  my $somatic_status = shift;

  my $row = qq{
        $header
        <td style="font-weight:bold">$source</td>
        <td>$version</td>
        <td style="max-width:800px">$desc</td>
        <td style="border-left:1px solid #CCC;padding:4px 2px">$data_type</td>
        <td style="text-align:center;width:22px;padding:2px 3px;border-left:1px solid #CCC">$phenotype</td>
        <td style="text-align:center;width:22px;padding:2px 3px;border-left:1px solid #DDD">$somatic_status</td>
  };
  return $row;
}


sub get_species_count {
  my $data_types = shift;
  my $species    = shift;
  my $database   = shift;
  my $hostname   = shift;

  my %count_by_type;
  
  foreach my $type (keys(%$data_types)) {
    next if ($type eq 'variation_set');
    my $sql = $data_types->{$type}{'count_spe'};

    my $sth = get_connection_and_query($database, $hostname, $sql);

    if ($sth) {
      while (my ($source_id,$count) = $sth->fetchrow_array) {
        $count_by_type{$type}{$source_id} = get_count($count);
      }
      $sth->finish();
    }
  }
  return \%count_by_type;
}

sub get_phenotype_types {
  my $data_types = shift;
  my $species    = shift;
  my $database   = shift;
  my $hostname   = shift;
  my $type       = 'phenotype_feature';

  my %phe_types;

  my $sql = $data_types->{$type}{'types'};

  my $sth = get_connection_and_query($database, $hostname, $sql);

  if ($sth) {
    while (my ($source_id,$types) = $sth->fetchrow_array) {
      $types =~ s/StructuralVariation/Structural Variation/;
      $phe_types{$type}{$source_id} = $types;
    }
    $sth->finish();
  }
  return \%phe_types;
}

sub get_species_set_count {
  my $param     = shift;
  my $species   = shift;
  my $database  = shift;
  my $hostname  = shift;
  my $data_type = 'variation_set';

  return '' if (!$data_type_example{$data_type});

  my $sql = $data_type_example{$data_type}{'count'};
  
  my $sth = get_connection_and_query($database, $hostname, $sql, [$param]);

  if ($sth) {
    my @result = $sth->fetchrow_array;
    return get_count($result[0]);
  }
  return '-'; 
}

sub get_count {
  my $count = shift;
  my $symbol = '+';
  
  my $count_label;
  my $count_display;
  my $bg_color;
  # From 1 to 9.9 million
  if ($count =~ /^(\d)(\d)\d{5}$/) {
    my $number = ($2!=0) ? "$1.$2" : $1;
    $count = "$number million";
    $count_label = "Over $count variants";
    $count_display = "$count$symbol";
    $bg_color = $colours{'few_million'};
  }
  # From 10 million
  elsif ($count =~ /^(\d+)\d{6}$/) {
    my $number = $1;
    $count = "$number million";
    $count_label = "Over $count variants";
    $count_display = "$count$symbol";
    $bg_color = $colours{'lot_million'};
  }
  # From 1,000 to 999,999
  elsif ($count =~ /^(\d+)\d{3}$/) {
    $count = "$1,000";
    $count_label = "Over $count variants";
    $count_display = "$count$symbol";
    $bg_color = $colours{'thousand'};
  }
  # From 1 to 999
  else {
    $count_label = "$count variants";
    $count_display = $count;
    $bg_color = $colours{'hundred'};
  }
  return qq{<span style="background-color:$bg_color;color:#FFF;border-radius:5px;padding:1px 2px;cursor:help;white-space:nowrap" title="$count_label"><small>$count_display</small></span>};
}

sub get_example {
  my $data_type  = shift;
  my $param      = shift;
  my $species    = shift;
  my $database   = shift;
  my $hostname   = shift;
  my $is_somatic = shift;

  return '' if (!$data_type_example{$data_type});

  my $sql_query = ($is_somatic) ? 'sql_som' : 'sql';
  my $sql = $data_type_example{$data_type}{$sql_query};
  my $url = $data_type_example{$data_type}{'url'};

  my $sth = get_connection_and_query($database, $hostname, $sql, [$param]);

  if ($sth) {
    my @result  = $sth->fetchrow_array;
    my $example = ($result[1] eq 'QTL') ? $result[2] : $result[0];
    my $url;
    if ($result[1] && $data_type_example{$data_type}{$result[1]}) {
      $url = $data_type_example{$data_type}{$result[1]};
    }
    else {
      $url = $data_type_example{$data_type}{'url'};
    }
    if ($example && $url) {
      $data_type =~ s/_/ /g;
      my $example_url = "/$species/$url$example";
      return qq{<a href="$example_url" target="_blank" title="See a $data_type example"><img src="$internal_link" alt="Link"/></a>};
    }
  }
  return '';
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
