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

# Script to generate an HTML page containing the phenotype variation sources of each species


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
my ($e_version,$html_file,$source_id,$source,$s_version,$s_description,$s_url,$s_type,$s_status,$s_data_types,$s_order,$hlist,$phost,$skip_name,$help);
my ($set_id,$set_name,$set_description);

## EG options
my ($site, $etype);

usage() if (!scalar(@ARGV));
 
GetOptions(
     'v=s'        => \$e_version,
     'o=s'        => \$html_file,
     'help!'      => \$help,
     'hlist=s'    => \$hlist,
     'phost=s'    => \$phost,
     'site=s'     => \$site,
     'skip_name!' => \$skip_name, 
     'etype=s'    => \$etype
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

my %colours = ( 'version' => '#090',
                'source'  => '#00F'
              );

my %colour_class = ( 'version'     => 'vdoc_new_version',
                     'source'      => 'vdoc_new_source',
                     'lot_million' => 'vdoc_million_1',
                     'few_million' => 'vdoc_million_2',
                     'thousand'    => 'vdoc_thousand', 
                     'hundred'     => 'vdoc_hundred'
                   );

my $phen_icon = '/i/val/var_phenotype_data_small.png';
my $internal_link = '/i/16/internal_link.png';

my %data_type_example = ('sql'       => qq{SELECT object_id, type, phenotype_id FROM phenotype_feature WHERE source_id=? AND type=? AND is_significant=1 AND type!="SupportingStructuralVariation" LIMIT 1},
                         'count'     => qq{SELECT COUNT(phenotype_id) FROM phenotype},
                         'count_spe' => qq{SELECT source_id, type, COUNT(phenotype_feature_id)
                                           FROM phenotype_feature WHERE type!="SupportingStructuralVariation" 
                                           GROUP BY source_id, type ORDER BY type ASC},
                         'Variation'           => 'Variation/Phenotype?v=',
                         'StructuralVariation' => 'StructuralVariation/Phenotype?v=',
                         'Gene'                => 'Gene/Phenotype?g=',
                         'QTL'                 => 'Phenotype/Locations?ph='                       
);

##############
### Header ###
##############
my $html_header = q{
<html>
<head>
  <title>Phenotype Sources</title>
  <script type="text/javascript">
    window.onload = function() {
      $('._ht').helptip({'track': true});
    };
  </script>
</head>

<body>

<div>
};


#############
### Title ###
#############
my $html_title = qq{
  <div style="float:left;width:75%">
    <h1 style="margin-top:15px">Ensembl Variation - Phenotype Sources Documentation</h1>

    <h2>List of sources providing phenotype/disease/trait associations for each species - $ecaption $e_version</h2>

    <div style="margin-bottom:20px">
      <a href="sources_documentation.html">See documentation for all the variant sources &rarr;</a>
    </div>
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
my %display_list;
my %top_display;

my $has_new_version = 0;
my $has_new_source  = 0;

my @new_status = ('version','source');

my $cols_sql2  = 'source_id, name, version, description, url, type, somatic_status, data_types';
my $term_sql2  = 'dbSNP';
my $condition1 = qq{SELECT $cols_sql2, 1 AS ordering FROM source WHERE name = "$term_sql2"};
my $condition2 = qq{SELECT $cols_sql2, 2 AS ordering FROM source WHERE name like "%$term_sql2%" AND name != "$term_sql2"};
my $condition3 = qq{SELECT $cols_sql2, 3 AS ordering FROM source WHERE name not like "%$term_sql2%"};

my $sql2 = qq{$condition1 UNION $condition2 UNION $condition3 ORDER BY ordering, name};
my $sql4 = qq{SELECT name, version FROM source};
my $sql_core = qq{SELECT meta_value FROM meta WHERE meta_key="species.display_name" LIMIT 1};


# Get the list of species and their common names
print STDERR "# Databases list:\n";
foreach my $hostname (@hostnames) {

  my $sql = qq{SHOW DATABASES LIKE '%variation_$e_version%'};
  my $sth = get_connection_and_query($database, $hostname, $sql);
  my $db_found = 0;
  
  # loop over databases
  while (my ($dbname) = $sth->fetchrow_array) {
    next if ($dbname =~ /^master_schema/ || $dbname =~ /^homo_sapiens_variation_\d+_37$/ || $dbname =~ /private/ || $dbname =~ /_variation_\d+_\d+_\w+$/ );
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
    
    # Get species display name
    my $display_name = $s_name;
    
    if (!$skip_name) {
      my $core_dbname = $dbname;
         $core_dbname =~ s/variation/core/i;
      my $sth_core = get_connection_and_query($core_dbname, $hostname, $sql_core);
      $display_name = $sth_core->fetchrow_array;  
      $display_name =~ s/saccharomyces/S\./i;
    }
    
    if ($top_species{$s_name}) {
      $top_display{$display_name} = 1;
    }
    $display_list{$display_name} = {'name' => $s_name, 'dbname' => $dbname, 'hostname' => $hostname};
  }
  die ("No variation databases found on $hostname for the version $e_version!") if ($db_found == 0);
}


print STDERR "\n# Get sources for each species:\n";
foreach my $display_name (sort { $top_display{$a} cmp $top_display{$b} || $a cmp $b} keys(%display_list)) {

  my $s_name   = $display_list{$display_name}{'name'};
  my $dbname   = $display_list{$display_name}{'dbname'};
  my $hostname = $display_list{$display_name}{'hostname'};
  
  print STDERR "- $display_name ($s_name) ...";
    
  # Get list of sources from the new databases
  my $sth2 = get_connection_and_query($dbname, $hostname, $sql2);
  $sth2->bind_columns(\$source_id,\$source,\$s_version,\$s_description,\$s_url,\$s_type,\$s_status,\$s_data_types,\$s_order);
    
 
  # Previous database (and sources)
  my $p_version = $e_version-1;
  my $sql3 = qq{SHOW DATABASES LIKE '%$s_name\_variation_$p_version%'};
  my $sth3 = get_connection_and_query($database, $previous_host, $sql3);
  my $p_dbname = $sth3->fetchrow_array;
 
  my %p_list;
  my $is_new_species = 0;
  if ($p_dbname) {
    # Previous sources
    my $sth4 = get_connection_and_query($p_dbname, $previous_host, $sql4);
    while (my @p = $sth4->fetchrow_array) {
      $p_list{$p[0]} = $p[1];
    }
  }
  else {
    $is_new_species = 1;
  }
    
  # Display the species at the top of the list
  if ($top_species{$s_name}) {
    $html_top_content .= source_phen_table($s_name,$display_name,$sth2,$is_new_species,$dbname,$hostname,\%p_list);
  }
  else {
    $html_content .= source_phen_table($s_name,$display_name,$sth2,$is_new_species,$dbname,$hostname,\%p_list);
  }
  print STDERR " done\n";
}

my $html_menu = create_menu();

if ($html_content ne '') {
  $html_top_content .= qq{
    <div style="background-color:#F0F0F0;margin:50px 0px 25px;padding:5px;border-top:2px solid #22949b;border-bottom:1px solid #22949b">
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

sub source_phen_table {
  my $name         = shift;
  my $display_name = shift;
  my $sth          = shift;
  my $is_new       = shift;
  my $db_name      = shift;
  my $hostname     = shift;
  my $p_list       = shift;

  
  my $s_name = ucfirst($name);
  my $species = $s_name;
     $species =~ s/_/ /;
  my $s_name_id = $name;
  
  my $count_phen = get_phenotype_count(\%data_type_example, $s_name, $db_name, $hostname);
  
  my $counts_species = get_phenotype_feature_count(\%data_type_example, $s_name, $db_name, $hostname);
  
  return '' if (!$counts_species || scalar(keys(%$counts_species))==0);
  
  if ($top_species{$name}) {
    $top_species_list{$name} = {name => $species, display => $display_name, s_name => $s_name, anchor => $s_name_id};
  }
  else {
    push (@species_list,{name => $species, display => $display_name, s_name => $s_name, anchor => $s_name_id});
  }


  my $phe_title = ($count_phen > 1) ? "See all ".lc($display_name)." phenotype entries" : "See ".lc($display_name)." phenotype entry";

  my $html = qq{<!-- $display_name ($species) -->};
  if ($is_new) {
    $html .= qq{
    <div style="padding-bottom:1px">
      <div style="float:left">
        <a href="/$s_name/Info/Index" title="$display_name Ensembl Home page" style="vertical-align:middle" target="_blank"><img src="/i/species/48/$s_name.png" alt="$display_name" class="sp-thumb" style="float:none;margin-right:0px;padding-right:0px;vertical-align:middle;background-color:#00F" /></a><h2 id="$s_name_id" style="display:inline;vertical-align:middle;margin-left:5px;padding:8px;background-color:#F0F0F0;color:#22949b">$display_name<span class="small vdoc_species_sci_name"> ($species)</span></h2><span style="padding-left:20px;color:#00F;font-weight:bold">New species!</span>
      </div>
      <div style="float:right;margin:25px 10px 0px 0px">
        <a href="/$s_name/Phenotype/All" title="$display_name Ensembl Phenotypes" style="vertical-align:middle" target="_blank"><img src="$phen_icon" style="border-radius:5px;border:1px solid #000;vertical-align:middle" alt="$phe_title" title="$phe_title" /></a>
        <span style="font-weight:bold;vertical-align:middle;margin-left:5px;color:#333" class="_ht ht" title="$count_phen phenotype(s)/disease(s)/trait(s) available for $display_name">$count_phen</span></span>
      </div>
      <div style="clear:both"></div>
    </div>
    };
  }
  else {
    $html .= qq{
    <div style="padding-bottom:3px">
      <div style="float:left">
        <a target="_blank" href="/$s_name/Info/Index" title="$display_name Ensembl Home page" style="vertical-align:middle"><img src="/i/species/48/$s_name.png" alt="$display_name" class="sp-thumb" style="float:none;margin-right:0px;padding-right:0px;vertical-align:middle;border-color:#22949b" /></a><h2 id="$s_name_id" style="display:inline;vertical-align:middle;margin-left:5px;padding:8px;background-color:#F0F0F0;color:#22949b">$display_name<span class="small vdoc_species_sci_name"> ($species)</span></h2>
      </div>
      <div style="float:right;margin:25px 10px 0px 0px">
        <a target="_blank" href="/$s_name/Phenotype/All" title="$display_name Ensembl Phenotypes" style="vertical-align:middle"><img src="$phen_icon" style="border-radius:5px;border:1px solid #000;vertical-align:middle" alt="$phe_title" title="$phe_title" /></a>
        <span style="font-weight:bold;vertical-align:middle;margin-left:5px;color:#333" class="_ht ht" title="$count_phen phenotype(s) available for $display_name">$count_phen</span></span>
      </div>
      <div style="clear:both"></div>
    </div>
    };
  }
  
  my $source_table;
  
  my $bg = 1;
  my @p_sources = keys(%{$p_list});
  my %other_flag;

  my $dt_class    = 'vdoc_data_types';
  my $type_class  = 'vdoc_type_style';
  my $count_class = 'vdoc_count_style';
  my $eg_class    = 'vdoc_example_style';
  my $spaces      = "        ";
  
  ########## Sources ##########    
  
  while ($sth->fetch) {
    
    # Check data types
    my %data_types = map { $_ => 1 } split(",", $s_data_types);
    next if (!$data_types{'phenotype_feature'});
  
    # Check if the source or its version is new
    my $s_new      = '';
    my $s_new_type = '';
    my $s_header   = '<td style="width:4px;padding:0px;margin:0px';
    
    # New source
    if (!grep {$_ eq $source} @p_sources) {
      $s_new_type = $new_status[1]; 
      $has_new_source = 1 if ($s_new_type);
    }
    # New version
    elsif ($p_list->{$source} ne $s_version){
      $s_new_type = $new_status[0]; 
      $has_new_version = 1 if ($s_new_type);
    }
   
    if ($s_new_type) {
      $species_news{$species}{$s_new_type} += 1;
      my $borders = ";border-top:1px solid #FFF;border-bottom:1px solid #FFF";
      $s_header .= $borders if ($bg == 1);
      $s_new = '<span style="color:'.$colours{$s_new_type}.'">New '.$s_new_type.'</span>' if ($s_new_type);
      $s_header .= ';background-color:'.$colours{$s_new_type};
    }
      
    $s_header .= '"></td>';

    # Source
    if ($s_url) {
      $source = qq{<a href="$s_url" style="text-decoration:none" target="_blank">$source</a>};
    }
    
    # Version
    $s_version = format_version($s_version);  
    
    # Somatic status
    my $is_somatic = ($s_status eq 'somatic') ? 1 : undef;
    
    # Data types
    my $data_type_string = '';
    my $counts;
    my $examples;

    foreach my $type (keys %{$counts_species->{$source_id}}) {
      my $type_label = ($type =~ /structural\s?variation/i) ? 'SV' : $type;
      my $dt_phe_title = "Provides $type phenotype association data";
      $data_type_string .= qq{\n$spaces<div class="$dt_class">\n$spaces  <div class="$type_class"><span class="_ht ht" title="$dt_phe_title">$type_label</span></div>};
 
      # Count
      my $count = $counts_species->{$source_id}{$type};
      $data_type_string .= qq{\n$spaces  <div class="$count_class">$count</div>};
      
      # Example     
      my $example = get_example($source_id, $s_name, $db_name, $hostname, $type);
      $data_type_string .= qq{\n$spaces  <div class="$eg_class">$example</div>};
      
      $data_type_string .= qq{\n$spaces  <div style="clear:both"></div>\n$spaces</div>};
    }
    $data_type_string = '-' if ($data_type_string eq '');
    
    $s_type = 'main' if (!defined($s_type));
    
    my $row = set_row($s_header,$source,$s_version,$s_description,$data_type_string);
    
    
    $source_table .= qq{
      <tr class="bg$bg">
        $row
      </tr>};
    if ($bg == 1) { $bg = 2; }
    else { $bg = 1; }
  }
 
  # Main source header
  $source_table = table_header('Source','main',\%other_flag).$source_table;
  
  $html .= qq{$source_table\n    </table>\n  </div>};
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
    $sth->execute(@$params);
  }
  else {
    $sth->execute;
  }
  return $sth;
}


sub create_menu {
  
  my %desc = ( 'version' => qq{New data version},
               'source'  => qq{New data source}
             );
  
  my $html = qq{
  <!-- Right hand side menu -->
  <div style="float:right">
    <div class="vdoc_menu vdoc_menu_phen">
      <div class="vdoc_menu_header vdoc_menu_header_phen">
        <img src="/i/16/info.png" style="vertical-align:top" alt="info" /> 
        Species with phenotype
      </div>
  };
  foreach my $top_sp (sort { $top_species{$a} <=> $top_species{$b} } keys(%top_species)) {
    my $species = $top_species_list{$top_sp};
    $html .= menu_list($species,\%desc) if (defined($species));
  }
  if (scalar(keys(%top_species_list))) {
    $html .= qq{<div class="vdoc_menu_separator"></div>};
  }
  
  foreach my $species (@species_list) {
    $html .= menu_list($species,\%desc);
  }
  
  my $v_colour  = $colour_class{'version'};
  my $s_colour  = $colour_class{'source'};
  my $lm_colour = $colour_class{'lot_million'};
  my $fm_colour = $colour_class{'few_million'};
  my $t_colour  = $colour_class{'thousand'};
  my $h_colour  = $colour_class{'hundred'};
  my $v_label   = $desc{'version'}.'(s)';
  my $s_label   = $desc{'source'}.'(s)';
  
  my $legend_div_id = 'legend';

  if ($has_new_version) {
    $html .= qq{    <div><span class="$v_colour vdoc_source_count" style="margin-left:5px">#</span><small> : $v_label</small></div>};
  }
  if ($has_new_source) {
    $html .= qq{    <div><span class="$s_colour vdoc_source_count" style="margin-left:5px">#</span><small> : $s_label</small></div>};
  }
  $html .= qq{
  </div>
  <div id="$legend_div_id" class="vdoc_menu vdoc_menu_phen">
    <!-- Legend header -->
    <div class="vdoc_menu_header vdoc_menu_header_phen">
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
        <td style="padding-top:6px;text-align:center;width:25px">
         <div class="vdoc_source_legend $v_colour"></div>
        </td>
        <td style="padding-top:4px"><b>New version</b> of the data<br />source in this release<br />for the species</td>
      </tr>
      <tr>
        <td style="padding-top:6px;text-align:center;width:25px">
          <div class="vdoc_source_legend $s_colour"></div>
        </td>
        <td style="padding-top:4px"><b>New data source</b> in this<br />release for the species</td>
      </tr>
    </table>
    
    <!-- Phenotype associations count colour legend -->
    <div style="border-top:1px dotted #22949b;margin-top:2px;padding:4px 0px 0px">
      <span style="padding-left:4px;font-weight:bold">Associations count:</span>
      <table>
        <tr>
          <td style="padding-top:4px;text-align:center">
            <span class="vdoc_count_legend $lm_colour"></span>
          </td>
          <td style="padding-top:4px">greater than 10 million</td>
        </tr>
        <tr>
          <td style="padding-top:4px;text-align:center">
            <span class="vdoc_count_legend $fm_colour"></span>
          </td>
          <td style="padding-top:4px">from 1 million to 9.9 million</td>
        </tr>
        <tr>
          <td style="padding-top:3px;text-align:center">
            <span class="vdoc_count_legend $t_colour"></span>
          </td>
          <td style="padding-top:3px">from 1,000 to 999,999</td>
        </tr>
        <tr>
          <td style="padding-top:3px;text-align:center">
            <span class="vdoc_count_legend $h_colour"></span>
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
          legend_element.style.right = "16px"; // Padding right of the page content
        }
        else {
          legend_element.style.position = "relative";
          legend_element.style.top = "auto";
          legend_element.style.right = "0px";
        }
      };
    </script>

  </div>
</div>
  };
  return $html;
}


sub menu_list {
  my $species = shift;
  my $desc = shift;

  my $name    = $species->{name};
  my $display = $species->{display};
  my $s_name  = $species->{s_name};
  my $anchor  = $species->{anchor};
  my $new_data = '';
  if ($species_news{$species->{name}}) {
    foreach my $type (@new_status) {
      next if (!$species_news{$species->{name}}{$type});
      my $count_type = $species_news{$species->{name}}{$type};
      my $label_colour = $colour_class{$type};
      my $label_desc = "$count_type ".lc($desc->{$type});
         $label_desc .= 's' if ($count_type > 1);

      my $label_pos = 'vdoc_new_'.$type.'_count';
      $new_data .= qq{<span class="_ht vdoc_source_count $label_colour $label_pos" title="$label_desc">$count_type</span>};
    }
  }
  my $img = $name;
  return qq{
  <div style="margin:0px 4px 5px">
    <div class="vdoc_menu_left">
      <img src="/i/species/16/$s_name.png" alt="$display" style="margin-right:4px;vertical-align:middle" />
      <a href="#$anchor" style="margin-right:3px;text-decoration:none;vertical-align:middle" title="$name">$display</a>
    </div>
    <div class="vdoc_menu_right">  
      $new_data
    </div>
    <div style="clear:both"></div>  
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


sub table_header {
  my $name = shift;
  my $type = shift;
  my $flag = shift;
  
  my $top_margin = ($type eq 'main') ? '6px' : '0px';

  my $border_color = qq{style="border-color:#22949b"};

  my $data_type_header = qq{
     <th class="vdoc_data_types vdoc_extra_column" style="background-color:#BBB">Phenotype/Disease/Trait
       <div class="vdoc_data_types">
         <div class="vdoc_type_style" style="text-align:center"><span class="_ht ht" $border_color title="Feature type association"><small>Type</small></span></div>
         <div class="vdoc_count_style" style="text-align:center"><span class="_ht ht" $border_color title="Phenotype associations count"><small>Count</small></span></div>
         <div class="vdoc_example_style"><span class="_ht ht" $border_color title="Example"><small>e.g.</small></span></div>
         <div style="clear:both"></div>
       </div>
     </th>
  };
  
  return qq{    <div style="margin:$top_margin 0px 30px">
      <table class="ss">
        <tr><th colspan="2">$name</th><th>Version</th><th style="max-width:800px">Description</th>$data_type_header</th></tr>
    };
}


sub set_row {
  my $header    = shift;
  my $source    = shift;
  my $version   = shift;
  my $desc      = shift;
  my $data_type = shift;
  my $phenotype = shift;

  my $row = qq{
        $header
        <td style="font-weight:bold">$source</td>
        <td>$version</td>
        <td style="max-width:800px">$desc</td>
        <td class="vdoc_extra_column">$data_type</td>
  };
  return $row;
}


sub get_phenotype_feature_count {
  my $data_types = shift;
  my $species    = shift;
  my $database   = shift;
  my $hostname   = shift;

  my %count_by_type;
  
  my $sql = $data_types->{'count_spe'};

  my $sth = get_connection_and_query($database, $hostname, $sql);

  if ($sth) {
    while (my ($source_id,$type,$count) = $sth->fetchrow_array) {
      $type =~ s/StructuralVariation/Structural Variation/;
      $count_by_type{$source_id}{$type} = get_count($count,$type);
    }
    $sth->finish();
  }
  return \%count_by_type;
}

sub get_phenotype_count {
  my $data_types = shift;
  my $species    = shift;
  my $database   = shift;
  my $hostname   = shift;

  my $count = 0;
  
  my $sql = $data_types->{'count'};

  my $sth = get_connection_and_query($database, $hostname, $sql);

  if ($sth) {
    $count = ($sth->fetchrow_array)[0];
    $sth->finish();
  }
  
  return thousandify($count);
}


sub get_count {
  my $count = shift;
  my $type  = lc(shift);
  my $symbol = '+';
  
  my $count_label;
  my $count_display;
  my $end_label = 'phenotype associations';
  my $bg_class;
  # From 1 to 9.9 million
  if ($count =~ /^(\d)(\d)\d{5}$/) {
    my $number = ($2!=0) ? "$1.$2" : $1;
    $count = "$number million";
    $count_label = "Over $count $type $end_label";
    $count_display = "$count$symbol";
    $bg_class = $colour_class{'few_million'};
  }
  # From 10 million
  elsif ($count =~ /^(\d+)\d{6}$/) {
    my $number = $1;
    $count = "$number million";
    $count_label = "Over $count $type $end_label";
    $count_display = "$count$symbol";
    $bg_class = $colour_class{'lot_million'};
  }
  # From 1,000 to 999,999
  elsif ($count =~ /^(\d+)\d{3}$/) {
    $count = "$1,000";
    $count_label = "Over $count $type $end_label";
    $count_display = "$count$symbol";
    $bg_class = $colour_class{'thousand'};
  }
  # From 1 to 999
  else {
    my $s = ($count > 1) ? 's' : '';
    $count_label = "$count $type phenotype association$s";
    $count_display = $count;
    $bg_class = $colour_class{'hundred'};
  }
  return qq{<span class="vdoc_dtype_count $bg_class" title="$count_label">$count_display</span>};
}

sub get_example {
  my $source_id  = shift;
  my $species    = shift;
  my $database   = shift;
  my $hostname   = shift;
  my $type       = shift;

  my $sql = $data_type_example{'sql'};
  my $url = $data_type_example{'url'};
  
  my $type_sql = $type;
     $type_sql =~ s/ //g;

  my @params = ($source_id,$type_sql);

  my $sth = get_connection_and_query($database, $hostname, $sql, \@params);
   
  if ($sth) {
    my @result  = $sth->fetchrow_array;
    my $example = ($result[1] eq 'QTL') ? $result[2] : $result[0];
    my $url;
    if ($result[1] && $data_type_example{$result[1]}) {
      $url = $data_type_example{$result[1]};
    }
    else {
      $url = $data_type_example{'url'};
    }
    if ($example && $url) {
      my $example_url = "/$species/$url$example";
      return qq{<a href="$example_url" target="_blank" title="See a $type example"><img src="$internal_link" alt="Link"/></a>};
    }
  }
  return '';
}

sub thousandify {
  my $value = shift;
  local $_ = reverse $value;
  s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
  return scalar reverse $_;
}

sub usage {
  
  print qq{
  Usage: perl phensources2html.pl [OPTION]
  
  Put all phenotype sources, for each species, into an HTML document.
  
  Options:

    -help           Print this message
      
    -v              Ensembl version, e.g. 65 (Required)
    -o              An HTML output file name (Required)      
    -phost          Host name (with port) where the previous databases are stored, e.g. ensembldb.ensembl.org:3306  (Required)
    -hlist          The list of host names (with port) where the new databases are stored, separated by a coma,
                    e.g. ensembldb.ensembl.org1:1234, ensembldb.ensembl.org2:1234 (Required)
    -skip_name      Flag to avoid the connection to the Core databases (use to retrieve the species display name, e.g. Human).
                    If the flag is used, the species name displayed will be the scientific name (e.g. Homo_sapiens).
    -site           The URL of the website (optional)
    -etype          The type of Ensembl, e.g. Plant (optional)
  } . "\n";
  exit(0);
}
