# Script to generate an HTML page containing the variation sources of each species

=head1 LICENSE

  Copyright (c) 1999-2013 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/legal/code_licence.html

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
my ($e_version,$html_file,$source_id,$source,$s_version,$s_description,$s_url,$s_type,$s_status,$hlist,$phost,$help);
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
if ($etype) {
  $ecaption .= ' '.ucfirst($etype);
}
# Settings
my $database = "";
my $login = "ensro";
my $pswd = "";
my $sep = "\t";
my $start = 0;
my %colours = ( 'version' => '#090', 'source'  => '#00F' );


##############
### Header ###
##############
my $html_header = qq{
<html>
<head>
  <title>Variation Sources</title>
</head>

<body>

<div>
};

my $html_title = qq{
  <h1 style="margin-top:15px">Ensembl Variation - Sources Documentation</h1>

  <h2>List of Variation sources for each species - $ecaption $e_version</h2>
};


##############
### Footer  ##
##############
my $html_footer = qq{
</div>
</body>
</html>};


my $html_content = '';
my @species_list;
my %species_news;

foreach my $hostname (@hostnames) {

  my $sql = qq{SHOW DATABASES LIKE '%variation_$e_version%'};
  my $sth = get_connection_and_query($database, $hostname, $sql);

  # loop over databases
  while (my ($dbname) = $sth->fetchrow_array) {
    next if ($dbname =~ /^master_schema/);
    print $dbname;
    $dbname =~ /^(.+)_variation/;
    my $s_name = $1;

    my %phe_sources;

    if ($etype) { # EG site - need to filter out species
      my $img_thumb = sprintf qq{eg-plugins/%s/htdocs/img/species/thumb_%s.png}, $etype, ucfirst($s_name);
      #  print "- checking for $img_thumb ... ";
      if (! -e $img_thumb) {
        print "\t... skipping \n";
        next;
      } 
    }
    print "\n";
    # Get list of sources from the new databases
    my $sql2 = qq{SELECT source_id, name, version, description, url, type, somatic_status FROM source};
    my $sth2 = get_connection_and_query($dbname, $hostname, $sql2);
    $sth2->bind_columns(\$source_id,\$source,\$s_version,\$s_description,\$s_url,\$s_type,\$s_status);
    
    
    # Previous database (and sources)
    my $p_version = $e_version-1;
    my $sql3 = qq{SHOW DATABASES LIKE '%$s_name\_variation_$p_version%'};
    my $sth3 = get_connection_and_query($database, $previous_host, $sql3);
    my $p_dbname = $sth3->fetchrow_array;
    
    my %p_list;
    my $is_new_species = 0;
    if ($p_dbname) {
      my $sql4 = qq{SELECT name, version FROM source};
      my $sth4 = get_connection_and_query($p_dbname, $previous_host, $sql4);
      while (my @p = $sth4->fetchrow_array) {
        $p_list{$p[0]} = $p[1];
      }
    }
    else {
      $is_new_species = 1;
    }
    
    # Phenotype
    my $sql_phe = qq{SELECT distinct source_id FROM phenotype_feature};
    my $sth_phe = get_connection_and_query($dbname, $hostname, $sql_phe);
    while (my @phe = $sth_phe->fetchrow_array) {
      $phe_sources{$phe[0]} = 1;
    }
    $sth_phe->finish();
    
    $html_content .= qq{\n  <br />\n} if ($start == 1);
    $html_content .= source_table($s_name,$sth2,$is_new_species,\%p_list,\%phe_sources);
    
    $start = 1 if ($start == 0);
  }
}

my $html_menu = create_menu();

## HTML/output file ##
open  HTML, "> $html_file" or die "Can't open $html_file : $!";
print HTML $html_header."\n";
print HTML $html_menu."\n";
print HTML $html_title."\n";
print HTML $html_content."\n";
print HTML $html_footer."\n";
close(HTML);


sub source_table {
  my $name        = shift;
  my $sth         = shift;
  my $is_new      = shift;
  my $p_list      = shift;
  my $phe_sources = shift;
  
  my $species = $name;
     $species =~ s/_/ /;
     $species =~ /^(\w)(.+)$/;
     $species = uc($1).$2;
  my $s_name = $species;
     $s_name =~ s/\s/_/g;
  my $s_name_id = lc($s_name);
  
  push (@species_list,{name => $species, anchor => $s_name_id});
  
  my $new_species = ($is_new == 1) ? qq{<span style="padding-left:20px;color:#00F;font-weight:bold">New species!</span>} : '';
  
  my $html = qq{
  <div id="$s_name_id" style="padding-left:0px;padding-bottom:3px">
    <img src="/i/species/48/$s_name.png" alt="$species" class="sp-thumb" style="float:none;margin-right:4px;vertical-align:middle" />
    <span style="font-weight:bold;font-size:1.1em;color:#333">$species</span>$new_species
  </div>
  };
  my $source_table;
  
  my $bg = 1;
  my @p_sources = keys(%{$p_list});
  my %other_flag;
  
  # Chip headers
  my $cbg = 1;
  my $chip_table;
     
  # LSDB headers
  my $lbg = 1;
  my $lsdb_table;
     
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

    
    # Display
    if ($s_url) {
      $source = qq{<a href="$s_url">$source</a>};
    }
    
    $s_version = format_version($s_version);
    
    # Somatic status
    my $s_somatic_status = somatic_status($s_status);
    
    # Phenotype
    my $phe_title = "Provides phenotype data";
    my $s_phenotype = ($phe_sources->{$source_id}) ? qq{<img src="phenotype.png" style="border-radius:5px;border:1px solid #000" alt="$phe_title" title="$phe_title" />} : '';
    
    # New source/version
    my $s_new_stuff = ($s_new_type) ? new_source_or_version($s_new_type) : '';
    
    $s_type = 'main' if (!defined($s_type));
    $other_flag{$s_type} = 1 if ($s_phenotype ne '' || $s_new_stuff ne '' || $s_somatic_status ne '-');
      
    my $new = ($s_new_stuff eq '') ? qq{<td style="border-left:1px solid #BBB;padding:0px"></td>} : qq{<td style="text-align:center;width:22px;padding:2px 3px;border-left:1px solid #BBB">$s_new_stuff</td>};    
    my $row = qq{
      $s_header
      <td>$source</td>
      <td>$s_version</td>
      <td>$s_description</td>
      $new
      <td style="text-align:center;width:22px;padding:2px 3px">$s_phenotype</td>
      <td style="text-align:center;width:22px;padding:2px 3px">$s_somatic_status</td>
    };
    #$row .= ($s_new_stuff eq '') ? qq{<td></td>} : qq{<td style="text-align:center;width:22px">$s_new_stuff</td>};
    
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
 
  # Main source header
  $source_table = table_header('Source','main',\%other_flag).$source_table;
  
  # Chip header
  $chip_table = table_header('Chip Source','chip',\%other_flag).$chip_table if ($chip_table);
     
  # LSDB header
  $lsdb_table = table_header('LSDB Source','lsdb',\%other_flag).$lsdb_table if ($lsdb_table);
  
  $html .= qq{$source_table\n    </table>\n  </div>};
  $html .= qq{$chip_table\n    </table>\n  </div>} if ($chip_table);
  $html .= qq{$lsdb_table\n    </table>\n  </div>} if ($lsdb_table);
  
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
  
  my ($host, $port) = split /\:/, $hname;

  # DBI connection 
  my $dsn = "DBI:mysql:$dbname:$host:$port";
  my $dbh = DBI->connect($dsn, $login, $pswd) or die "Connection failed";

  my $sth = $dbh->prepare($sql);
  $sth->execute;
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
  <div style="margin-left:8px;margin-top:2px;padding-bottom:2px;background-color:#F2F2F2;color:#333;border-radius:5px">
    <div style="padding:5px;font-weight:bold;color:#FFF;background-color:#336;border-top-left-radius:5px;border-top-right-radius:5px">
      <img src="/i/16/rev/info.png" style="vertical-align:top" />
      Species list
    </div> 
    <ul style="padding:5px 5px 5px 25px;margin-bottom:0px">
  };
  foreach my $species (@species_list) {
    my $name = $species->{name};
    my $anchor = $species->{anchor};
    my $new_data = '';
    if ($species_news{$species->{name}}) {
      my @types = sort {$b cmp $a} keys(%{$species_news{$species->{name}}});
      foreach my $type (@types) {
        my $label_colour = $colours{$type};
        my $label_desc = $desc{$type};
        $new_data .= qq{<span style="$label_style;margin-right:5px;background-color:$label_colour" title="$label_desc"></span>};
      }
    }
    #my $new_data = ($species_news{$species->{name}}) ? qq{<span style="margin-left:10px;display:inline-block;height:10px;width:10px;border-radius:5px;background-color:#0A0"></span>} : '';
    $html .= qq{\n      <li><a href="#$anchor" style="margin-right:5px">$name</a>$new_data</li>};
  }
  my $v_colour = $colours{'version'};
  my $s_colour = $colours{'source'};
  my $v_label  = $desc{'version'};
  my $s_label  = $desc{'source'};
  
  $html .= sprintf ( qq{
    </ul>
    <span style="$label_style;margin-left:5px;background-color:$v_colour"></span><small> : $v_label</small>
    <br />
    <span style="$label_style;margin-left:5px;background-color:$s_colour"></span><small> : $s_label</small>
  </div>
  <br />
  <div id="legend" style="margin-left:8px;margin-top:2px;background-color:#F2F2F2;color:#333;border-radius:5px">
    <div style="padding:5px;font-weight:bold;color:#FFF;background-color:#336;border-top-left-radius:5px;border-top-right-radius:5px">
      <img src="/i/16/rev/info.png" style="vertical-align:top" />
      Icons legend
    </div> 
    <table>
      <tr>
        <td style="padding-top:8px">%s</td>
        <td style="padding-top:6px"><b>New version</b> of the data source<br />in this release for the species</td>
      </tr>
      <tr>
        <td style="padding-top:8px">%s</td>
        <td style="padding-top:6px"><b>New data source</b> in this<br />release for the species</td>
      </tr>
      <tr>
        <td style="padding-top:6px"><img src="phenotype.png" style="border-radius:5px;border:1px solid #000;margin-right:1px" alt="Provides phenotype data" title="Provides phenotype data"/></td>
        <td style="padding-top:6px">Source which provides<br />phenotype association data</td>
      </tr>
      <tr>
        <td style="padding-top:6px">%s</td>
        <td style="padding-top:6px">The source contains only<br />germline data</td>
      </tr>
      <tr>
        <td style="padding-top:6px">%s</td>
        <td style="padding-top:6px">The source contains only<br />somatic data</td>
      </tr>
      <tr>
        <td style="padding-top:6px">%s</td>
        <td style="padding-top:6px">The source contains both<br />germline and somatic data</td>
      </tr>
    </table>
  </div>
</div>
  },
  new_source_or_version('version'),
  new_source_or_version('source'),
  somatic_status('germline'),
  somatic_status('somatic'),
  somatic_status('mixed'));
  return $html;
}

sub new_source_or_version {
  my $type = shift;
  my $label = ($type eq 'version') ? 'V' : 'S';
  my $color = $colours{$type};
  
  return qq{
       <div style="color:$color;font-size:0.8em;text-align:center;margin:0px;padding:0px"><span style="margin:0px;padding:0px">New</span> <span style="margin:0px;padding:0px">$type</span></div>
     };
}

sub somatic_status {
  my $type = shift;
  my $html;
  if ($type eq 'germline') {
     $html .= qq{
       <div style="border-radius:5px;border:1px solid #000;width:20px;height:20px;background-color:#00C;" title="$type data"></div>
     };
  }
  elsif ($type eq 'somatic') {
    $html .= qq{
      <div style="border-radius:5px;border:1px solid #000;width:20px;height:20px;background-color:#C00;" title="$type data"></div>
    };
  }
  elsif ($type eq 'mixed') {
    $html .= qq{
      <div style="border-radius:5px;border:1px solid #000;width:20px;height:20px;background-color:#00C;" title="$type data">
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
  
  my $header_col;
  if ($flag->{$type}) {
    my $alt_text = qq{Click here to go to the icons table description};
    $header_col = qq{
    <th colspan=3 style="width:56px;text-align:center;border-left:1px solid #BBB;background-color:#BBB">
      Other<a href="#legend" class="ht _ht" title="$alt_text"><img src="/i/16/info.png" style="position:relative;top:2px;width:12px;height:12px;margin-left:3px" alt="$alt_text"/></a>
    </th>}; 
  }
  else {  
    $header_col = qq{<th colspan=3></th>};
  }
  
  return qq{
    <br />
    <div>
      <table class="ss" style="width:75%">
        <tr><th colspan="2">$name</th><th>Version</th><th>Description</th></th>$header_col</tr>
    };
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
