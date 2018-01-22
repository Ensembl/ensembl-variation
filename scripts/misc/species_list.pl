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
my ($e_version,$html_file,$hlist,$phost,$user,$help);
## EG options
my ($site, $etype);

usage() if (!scalar(@ARGV));
 
GetOptions(
     'v=s'     => \$e_version,
     'o=s'     => \$html_file,
     'help!'   => \$help,
     'hlist=s' => \$hlist,
     'phost=s' => \$phost,
     'user=s'  => \$user,
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
my $p_version = $e_version-1;
my $detailed_counts = 'species_detailed_counts.html';

my $html;
   
my %colours = ( 'lot_million' => { 'order' => 1, 'colour' => 'vdoc_million_1', 'legend' => 'From 10 million'},
                'few_million' => { 'order' => 2, 'colour' => 'vdoc_million_2', 'legend' => 'From 1 million to 9.9 million'},
                'thousand'    => { 'order' => 3, 'colour' => 'vdoc_thousand',  'legend' => 'From 1,000 to 999,999'},
                'hundred'     => { 'order' => 4, 'colour' => 'vdoc_hundred',   'legend' => 'From 1 to 999'}
              );              
              
my %tables = ( 'Sample'             => { 'order' => 2 , 'anchor' => '#genotype',             'table' => 'compressed_genotype_var'},
               'Population'         => { 'order' => 3 , 'anchor' => '#genotype',             'table' => 'population_genotype'},
               'Phenotype'          => { 'order' => 4 , 'anchor' => '#phenotype',            'table' => 'phenotype_feature'},
               'Citation'           => { 'order' => 5 , 'anchor' => '#citation',             'table' => 'variation_citation'},
               'Structural variant' => { 'order' => 1 , 'anchor' => '#structural_variation', 'table' => 'structural_variation'}
             );
             
my %columns = ( 'SIFT'     => {'order' => 1 , 'anchor' => '#prediction', 'table' => 'meta', 'column' => 'meta_key', 'value' => 'sift_version'},
                'PolyPhen' => {'order' => 2 , 'anchor' => '#prediction', 'table' => 'meta', 'column' => 'meta_key', 'value' => 'polyphen_version'}
              );
                        
my %species_list;
my %display_list;

my $sql   = qq{SHOW DATABASES LIKE '%$db_type\_$e_version%'};
my $sql2  = qq{SELECT COUNT(variation_id) FROM variation};
my $sql2a = qq{SELECT DISTINCT s.name FROM variation v, source s WHERE s.source_id=v.source_id ORDER BY s.name};
my $sql_core = qq{SELECT meta_value FROM meta WHERE meta_key="species.display_name" LIMIT 1};

foreach my $hostname (@hostnames) {
  
  my $sth = get_connection_and_query($database, $hostname, $sql);

  # loop over databases
  while (my ($dbname) = $sth->fetchrow_array) {
    next if ($dbname =~ /^master_schema/ || $dbname =~ /^homo_sapiens_variation_\d+_37$/ || $dbname =~ /private/ || $dbname =~ /_variation_\d+_\d+_\w+$/ );
    
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
    
    # Get species display name
    my $core_dbname = $dbname;
       $core_dbname =~ s/variation/core/i;
    my $sth_core = get_connection_and_query($core_dbname, $hostname, $sql_core);
    my $display_name = $sth_core->fetchrow_array;  
       $display_name =~ s/saccharomyces/S\./i;
    $species_list{$s_name}{'name'} = $display_name;
    $display_list{$display_name} = $s_name;
    
    # Count the number of variations
    my $sth2 = get_connection_and_query($dbname, $hostname, $sql2);
    my $count_var = $sth2->fetchrow_array;
    $sth2->finish;
    $species_list{$s_name}{'count'} = round_count($count_var);

    # Count the number of variantion sources
    my $sth2a = get_connection_and_query($dbname, $hostname, $sql2a);
    $species_list{$s_name}{'sources'} = [];
    while(my $src = ($sth2a->fetchrow_array)[0]) {
      push(@{$species_list{$s_name}{'sources'}},$src);
    }
    $sth2a->finish;

    # Previous database
    my $sql3 = qq{SHOW DATABASES LIKE '%$s_name\_variation_$p_version%'};
    my $sth3 = get_connection_and_query($database, $phost, $sql3);
    my $p_dbname = $sth3->fetchrow_array;

    if ($p_dbname) {
      # Previous sources
      my $sth4 = get_connection_and_query($p_dbname, $phost, $sql2);
      my $count_p_var = $sth4->fetchrow_array;
      $sth4->finish;
      $species_list{$s_name}{'p_count'} = round_count_diff($count_var-$count_p_var);
    }
  }
}

my $count_species = scalar(keys(%species_list));

# Get the populated tables by species
my $species_data_tables = get_species_data_tables();

my $species_data_columns = get_species_data_columns();

my $th_bg = qq{background-color:#BBB};
my $th_border_left = qq{border-left:1px solid #DDD};
my $th_border_left_top = qq{style="$th_border_left;text-align:center"};




my $data_tables_header = '';
my $header_link = qq{<a class="_ht" style="text-decoration:none" title="See detailed counts" href="};
foreach my $column_title (sort { $tables{$a}{'order'} <=> $tables{$b}{'order'} } keys(%tables)) {
  my $link = $header_link.$detailed_counts.$tables{$column_title}{'anchor'}.qq{">$column_title</a>};
  $data_tables_header .= qq{<th style="$th_bg;$th_border_left">$link</th>};
}
my $data_columns_header = '</th>';
foreach my $column_title (sort { $columns{$a}{'order'} <=> $columns{$b}{'order'} } keys(%columns)) {
  my $link = $header_link.$detailed_counts.$columns{$column_title}{'anchor'}.qq{">$column_title</a>};
  $data_columns_header .= qq{<th style="$th_bg;$th_border_left">$link</th>};
}

my $html_content = qq{
  <table class="ss" style="width:auto">
    <tr class="ss_header">
      <th></th>
      <th $th_border_left_top colspan="3">Short variant</th>
      <th $th_border_left_top colspan="1">Long variant</th>
      <th $th_border_left_top colspan="2">Genotype</th>
      <th $th_border_left_top colspan="2">Association</th>
      <th $th_border_left_top colspan="2">Prediction</th>
    </tr>
    
    <tr class="ss_header">
      <th style="$th_bg">Species</th>  
      <th style="$th_bg;$th_border_left">Sequence variant</th>
      <th style="$th_bg;padding-left:0px">
        <span class="_ht ht" title="Sequence variant count difference with the previous Ensembl release (v.$p_version)">
          <small>(e!$p_version &rarr; e!$e_version)</small>
        </span>
      </th>
      <th style="$th_bg;$th_border_left">Source(s)</th>
      $data_tables_header
      $data_columns_header
    </tr>};
my $bg = '';

foreach my $display_name (sort keys(%display_list)) {

  my $sp = $display_list{$display_name};
  my $label = $species_list{$sp}{'label'};
  my $uc_sp = ucfirst($sp);      
  my $img_src = "/i/species/48/$uc_sp.png";
  my $display_name = $species_list{$sp}{'name'};
  my $var_count = $species_list{$sp}{'count'};
  my $var_p_count = $species_list{$sp}{'p_count'};
  
  # Source(s)
  my @sources_list = sort { $a !~ /^dbSNP$/ cmp $b !~ /^dbSNP$/ || $a cmp $b } @{$species_list{$sp}{'sources'}};
  my $var_src_count  = @sources_list;
     $var_src_count .= ($var_src_count == 1) ? ' source' : ' sources';
  my $var_src_title  = ($var_src_count == 1) ? $sources_list[0] : join(' , ', @sources_list);
  
  $html_content .= qq{
  <tr$bg style="vertical-align:middle">
    <td>
      <div>
        <div style="float:left">
          <a href="/$uc_sp/Info/Index" title="$label Ensembl Home page" style="vertical-align:middle" target="_blank">
            <img src="$img_src" alt="$label" class="sp-thumb" style="vertical-align:middle;width:28px;height:28px" />
          </a>
         </div>
         <div style="float:left;margin-left:4px">
           <div class="bigtext">$display_name</div>
           <div class="small" style="font-style:italic">$label</div>
         </div>
         <div style="clear:both"></div>
       </div>
    </td>
    <td style="text-align:right">$var_count</td>
    <td style="text-align:right">$var_p_count</td>
    <td style="text-align:right">
      <a href="sources_documentation.html#$sp" class="ht _ht" style="text-decoration:none" title="$var_src_title">$var_src_count</a>
    </td>\n};
  
  
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


# Legend
my $html_legend = qq{
<span style="border:1px #DDD solid;padding:4px">
  <span style="margin-right:5px;font-weight:bold">Colour legend: </span>
};
foreach my $type (sort { $colours{$a}{'order'} <=> $colours{$b}{'order'} } keys(%colours)) {
  my $desc  = $colours{$type}{'legend'};
  my $class = $colours{$type}{'colour'};
  $html_legend .= qq{  
  <span style="margin-left:20px">
    <span class="vdoc_count_legend $class" style="margin-right:5px"></span>
    <span>$desc</span>
  </span>};
}
$html_legend .= qq{
</span>
};


## HTML/output file ##
open  HTML, "> $html_file" or die "Can't open $html_file : $!";
print HTML qq{<p style="padding-top:0px;margin-top:0px">There are currently <span style="font-weight:bold;font-size:1.1em;color:#000">$count_species</span> variation databases in Ensembl:</p>\n};
print HTML $html_content;
print HTML $html_legend;
print HTML qq{<p style="padding-top:15px">The <b>full list of species</b> with their assembly versions in Ensembl is available <a href="/info/about/species.html">here</a>.</p>\n};
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
  my $count = shift;
  my $type = 'variants';
  my $symbol = '+';
  
  my $count_label;
  my $count_display;
  my $bg_class;
  # From 1 to 9.9 million
  if ($count =~ /^(\d)(\d)\d{5}$/) {
    my $number = ($2!=0) ? "$1.$2" : $1;
    $count = "$number million";
    $count_label = "Over $count $type";
    $count_display = "$count$symbol";
    $bg_class = $colours{'few_million'}{'colour'};
  }
  # From 10 million
  elsif ($count =~ /^(\d+)\d{6}$/) {
    my $number = $1;
    $count = "$number million";
    $count_label = "Over $count $type";
    $count_display = "$count$symbol";
    $bg_class = $colours{'lot_million'}{'colour'};
  }
  # From 1,000 to 999,999
  elsif ($count =~ /^(\d+)\d{3}$/) {
    $count = "$1,000";
    $count_label = "Over $count $type";
    $count_display = "$count$symbol";
    $bg_class = $colours{'thousand'}{'colour'};
  }
  # From 1 to 999
  else {
    $count_label = "$count $type";
    $count_display = "$count";
    $bg_class = $colours{'hundred'}{'colour'};
  }
  return qq{<span class="vdoc_var_count $bg_class" title="$count_label">$count_display</span>};
}

sub round_count_diff {
  my $count = shift;
  my $type = 'variants';

  my ($count_label,$colour,$symbol,$label);

  if ($count == 0) {
    return '-';
  }
  elsif ($count > 0) {
    $colour = ' style="color:#090"';
    $symbol = '+';
    $label  = 'more'
  }
  else {
    $colour = '';
    $symbol = '-';
    $label  = 'less';
  }
  $count = abs($count);

  # From 1 to 9.9 million
  if ($count =~ /^(\d)(\d)\d{5}$/) {
    my $number = ($2!=0) ? "$1.$2" : $1;
    $count = "$number million";
    $count_label = "Over $count $label $type";
  }
  # From 10 million
  elsif ($count =~ /^(\d+)\d{6}$/) {
    my $number = $1;
    $count = "$number million";
    $count_label = "Over $count $label $type";
  }
  # From 1,000 to 999,999
  elsif ($count =~ /^(\d+)\d{3}$/) {
    $count = "$1,000";
    $count_label = "Over $count $label $type";
  }
  # From 1 to 999
  else {
    $count =~ /(\d+)$/;
    $count = $1;
    $count_label = "$count $label $type";
  }
  return qq{<span$colour title="$count_label"><small>($symbol$count)</small></span>};
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
  Usage: perl species_list.pl [OPTION]
  
  Put data information for each species, into an HTML document.
  
  Options:

    -help           Print this message
      
    -v              Ensembl version, e.g. 65 (Required)
    -o              An HTML output file name (Required)      
    -hlist          The list of host names (with port) where the new databases are stored, separated by a coma,
                    e.g. ensembldb.ensembl.org1:1234, ensembldb.ensembl.org2:1234 (Required)
    -phost          Host name (with port) where the previous databases are stored, e.g. ensembldb.ensembl.org:3306  (Required)
    -user           MySQL user name (Required)
    -site           The URL of the website (optional)
    -etype          The type of Ensembl, e.g. Plant (optional)
  } . "\n";
  exit(0);
}
