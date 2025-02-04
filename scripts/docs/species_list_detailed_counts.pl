# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2025] EMBL-European Bioinformatics Institute
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
  <http://www.ensembl.org/Help/Contact>.

=cut


use Bio::EnsEMBL::Registry;
use DBI;
use strict;
use POSIX;
use Getopt::Long;
use JSON;
use File::Basename;

my $dirname = dirname(__FILE__);
require "$dirname/utils.pl";

###############################################################
##########             CONFIGURE                        #######
###############################################################
my ($e_version,$html_file,$hlist,$phost,$user,$port,$config,$d_dir,$p_data,$skip_prediction,$help);
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
     'port=i'  => \$port,
     'config=s' => \$config,
     'd_dir=s'  => \$d_dir,
     'p_data=s' => \$p_data,
     'site=s'  => \$site,
     'etype=s' => \$etype,
     'skip_prediction' => \$skip_prediction
);

if (!$e_version) {
  print STDERR "> Error! Please give an Ensembl version, using the option '-v' \n";
  usage();
}
if (!$html_file) {
  print STDERR "> Error! Please give an output file using the option '-o'\n";
  usage();
}
if (!$hlist) {
  print STDERR "> Error! Please give the list of host names where the new databases are stored using the option '-hlist'\n";
  usage();
}
if (!$user) {
  print STDERR "> Error! Please give user name using the option '-user'\n";
  usage();
}

usage() if ($help);

my $server_name = 'http://static.ensembl.org';
my $ecaption = 'Ensembl';
my @hostnames = split /,/, $hlist;

if ($site) {
  $server_name = $site;
}

# Get the vcf config file location
my $vcf_config_file = $dirname . '/../../modules/Bio/EnsEMBL/Variation/DBSQL/vcf_config.json';

if ($config){
  $vcf_config_file = $config;
}

# Get the local dir where the vcf files are located
my $data_dir = "/nfs/production/flicek/ensembl/production/ensemblftp/data_files/vertebrates";

if ($d_dir){
  $data_dir = $d_dir;
}

# Settings
my $database = "";
my $pswd = "";
my $db_type = 'variation';
my $default_port = 3306;
$port ||= $default_port;
my $p_version = $e_version-1;

my $html;

# Header
my $html_header = q{
<html>
<head>
  <title>Detailed species data count</title>
  <script type="text/javascript">
    window.onload = function() {
      $('.ht').helptip({'track': true});
    };
  </script>
</head>

<body>

<div>
};

# Title
my $html_title = qq{
  <div style="float:left;width:80%">
    <h1 style="margin-top:15px">Ensembl Variation - Detailed species data count</h1>

    <h2>List of data counts by category and species - $ecaption $e_version</h2>

    <div>
      <a href="sources_documentation.html">See documentation for all the variant sources &rarr;</a>
    </div>
};

# Legend
my $html_legend = qq{
  <!-- Right hand side legend -->
  <div style="float:right;max-width:220px;top:20px">
    <div class="vdoc_menu vdoc_menu_phen">
      <div class="vdoc_menu_header vdoc_menu_header_phen">
        <img src="/i/16/info.png" style="vertical-align:top" alt="info" /> 
        Colour legend
      </div>
      <table>
};

# Footer 
my $html_footer = q{
  </div>
</div>
</body>
</html>};

my $prediction = 'Prediction';
   
my %colours = ( 'hundred_million' => { 'order' => 1, 'colour' => 'vdoc_million_0', 'legend' => 'From 100 million'},
                'lot_million'     => { 'order' => 2, 'colour' => 'vdoc_million_1', 'legend' => 'From 10 to 99.9 million'},
                'few_million'     => { 'order' => 3, 'colour' => 'vdoc_million_2', 'legend' => 'From 1 million to 9.9 million'},
                'thousand'        => { 'order' => 4, 'colour' => 'vdoc_thousand',  'legend' => 'From 1,000 to 999,999'},
                'hundred'         => { 'order' => 5, 'colour' => 'vdoc_hundred',   'legend' => 'From 1 to 999'},
                'zero'            => { 'order' => 6, 'colour' => 'vdoc_zero',      'legend' => 'No data'}
              );  

my $sql = qq{SHOW DATABASES LIKE '%$db_type\_$e_version%'};

my @genotype_projects = ('1000 Genomes', 'gnomAD', 'TOPMed', 'UK10K', 'Mouse Genomes', 'NextGen');
my $genotypes_list = qq{<ul><li>}.join(' Project</li><li>',@genotype_projects).qq{ Project</li></ul>};

my %sql_list = ( "Structural variant" => { 'sqla'   => { 'sql'   => q{SELECT COUNT(sv.structural_variation_id) FROM structural_variation sv, source s 
                                                                      WHERE sv.is_evidence=0 AND s.source_id=sv.source_id AND s.name IN ("DGVa", "dbVar")},
                                                         'label' => 'Structural variant'
                                                       }, 
                                           'sqlb'   => { 'sql'   => q{SELECT COUNT(sv.structural_variation_id) FROM structural_variation sv, source s 
                                                                     WHERE sv.is_evidence=1 AND s.source_id=sv.source_id AND s.name IN ("DGVa", "dbVar")},
                                                         'label' => 'Supporting evidence'
                                                       },
                                         },
                   "Citation"         => { 'sqla'   => { 'sql'   => q{SELECT COUNT(DISTINCT variation_id) FROM variation_citation},
                                                         'label' => 'Cited sequence variant'
                                                       }, 
                                           'sqlb'   => { 'sql'   => q{SELECT COUNT(DISTINCT publication_id) FROM variation_citation},
                                                         'label' => 'Publication'
                                                       },
                                         },
                   "Phenotype"        => { 'sqla'   => { 'sql'   => q{SELECT COUNT(phenotype_id) FROM phenotype},
                                                         'label' => 'Phenotype'
                                                       },
                                           'sqlb'   => { 'sql'   => q{SELECT COUNT(phenotype_feature_id) FROM phenotype_feature},
                                                         'label' => 'Phenotype association'
                                                       },
                                           'extra'  => q{The list of phenotype/disease/trait association sources by species is available in the page <a href="../phenotype/sources_phenotype_documentation.html">Phenotype sources</a>.}
                                         },
                   "Genotype"         => { 'sqla'   => { 'sql'   => q{SELECT COUNT(distinct variation_id) FROM compressed_genotype_var},
                                                         'label' => 'Variants with sample genotype'
                                                       },
                                           'sqlb'   => { 'sql'   => q{SELECT COUNT(distinct variation_id) FROM population_genotype},
                                                         'label' => 'Variants with population genotype'
                                                       },
                                           'vcf'    => 1
                                         },
                   $prediction        => { 'sqla'  => { 'sql'   => q{SELECT COUNT(distinct vf.variation_id) FROM variation_feature vf, transcript_variation tv, meta m 
                                                                     WHERE vf.variation_feature_id=tv.variation_feature_id AND m.meta_key="sift_version" 
                                                                     AND tv.sift_score IS NOT NULL},
                                                        'label' => 'Variants with SIFT data'
                                                      },
                                           'sqlb' =>  { 'sql'   => q{SELECT COUNT(distinct vf.variation_id) FROM variation_feature vf, transcript_variation tv, meta m
                                                                     WHERE vf.variation_feature_id=tv.variation_feature_id AND m.meta_key="polyphen_version" 
                                                                     AND tv.polyphen_score IS NOT NULL},
                                                        'label' => 'Variants with PolyPhen data'
                                                      },
                                         } # Too long ?
               );
my @type_order = ("Structural variant","Genotype","Phenotype","Citation",$prediction);

my $sql_core = qq{SELECT meta_value FROM meta WHERE meta_key="species.display_name" LIMIT 1};

# Get the populated tables by species
my $bg = '';
my $th_bg = qq{background-color:#BBB};
my $th_border_left = qq{border-left:1px solid #DDD};
my $th_border_left_top = qq{style="$th_border_left;text-align:center"};

# read config from JSON config file
open IN, $vcf_config_file or throw("ERROR: Could not read from config file $vcf_config_file\n");
local $/ = undef;
my $json_string = <IN>;
close IN;

my $vcf_config = JSON->new->decode($json_string) or throw("ERROR: Failed to parse config file $vcf_config_file\n");

# read previous release data
my $prev_data;
if ($p_data) {
  open IN, $p_data or throw("ERROR: Could not read from config file $p_data");
  local $/ = undef;
  $json_string = <IN>;
  close IN;

  $prev_data = JSON->new->decode($json_string) or throw("ERROR: Failed to parse config file $p_data");
}

###############################################################
##########             MAIN PART                       ########
###############################################################

foreach my $type (@type_order) {

  if ($skip_prediction && $type eq $prediction) {
    print STDERR "'skip_prediction' option used: the category $prediction will be skipped.\n";
    next;
  }

  print STDERR "# $type\n";

  if (!$sql_list{$type}) {
    print STDERR "Can't recognise the category '$type'! Skip this category.\n";
    next;
  }

  print STDERR "\n# $type\n";
  my $lc_type = lc($type);
  
  my $anchor = $lc_type;
     $anchor =~ s/ /_/g;
     
  my %species_list;
  my %display_list;
  
  # Column A
  my $sql_a   = $sql_list{$type}{'sqla'}{'sql'};
  my $label_a = $sql_list{$type}{'sqla'}{'label'};
  
  # Column B
  my $sql_b   = $sql_list{$type}{'sqlb'}{'sql'};
  my $label_b = $sql_list{$type}{'sqlb'}{'label'};
  my $b_type;
  
  foreach my $hostname (@hostnames) {
    
    my $sth = get_connection_and_query($database, $hostname, $sql);

    # loop over databases
    while (my ($dbname) = $sth->fetchrow_array) {
      next if ($dbname !~ /^[a-z][a-z_]*_[a-z]+_variation_\d+_\d+$/i);
      next if ($dbname =~ /^(master_schema|drosophila|saccharomyces|ciona)/ || $dbname =~ /^homo_sapiens_variation_\d+_37$/ || $dbname =~ /private/);
      print STDERR "${dbname}\n";
      
      $dbname =~ /^(.+)_variation/;
      my $s_name = $1;
      
      if ($etype) { # EG site - need to filter out species
        my $img_thumb = sprintf qq{eg-plugins/%s/htdocs/img/species/thumb_%s.png}, $etype, ucfirst($s_name);
        #  print STDERR "- checking for $img_thumb ... ";
        if (! -e $img_thumb) {
          print STDERR "\t... skipping \n";
          next;
        } 
      }
      print STDERR "\n";
      
      my $label_name = ucfirst($s_name);
         $label_name =~ s/_/ /g;
      $species_list{$s_name}{label} = $label_name;
      
      # Get species display name
      my $core_dbname = $dbname;
         $core_dbname =~ s/variation/core/i;
      my $sth_core = get_connection_and_query($core_dbname, $hostname, $sql_core);
      my $display_name = $sth_core->fetchrow_array;
      $species_list{$s_name}{'name'} = $display_name;
      $display_list{$display_name} = $s_name;
      
      # Count the number of variations - column A
      my $sth_a = get_connection_and_query($dbname, $hostname, $sql_a);
      my $res_a = $sth_a->fetchrow_array;
      $sth_a->finish;
      $species_list{$s_name}{'a'} = $res_a;

      # Count the number of variations - column B
      my $sth_b = get_connection_and_query($dbname, $hostname, $sql_b);
      my $res_b = $sth_b->fetchrow_array;
      $sth_b->finish;
      $species_list{$s_name}{'b'} = $res_b;

      # Previous database
      my $sql3 = qq{SHOW DATABASES LIKE '%$s_name\_variation_$p_version%'};
      my $sth3 = get_connection_and_query($database, $phost, $sql3);
      my $p_dbname = $sth3->fetchrow_array;

      if ($p_dbname) {

        # Previous variants - column A
        my $sth4a = get_connection_and_query($p_dbname, $phost, $sql_a);
        my $p_res_a = $sth4a->fetchrow_array;
        $sth4a->finish;
        $species_list{$s_name}{'p_a'} = $res_a-$p_res_a;

        # Previous variants - column B
        my $p_res_b;
        if ($res_b =~ /^\d+$/) {
          $b_type = "num";
          my $sth4b = get_connection_and_query($p_dbname, $phost, $sql_b);
          $p_res_b = $sth4b->fetchrow_array;
          $sth4b->finish;
          $species_list{$s_name}{'p_b'} = $res_b-$p_res_b;
        }
      }
      
      # Add vcf data if exist
      if ($sql_list{$type}{'vcf'}){
        foreach my $project (@{ $vcf_config->{'collections'} }) {
          next if $project->{annotation_type} eq 'cadd' || $project->{annotation_type} eq 'gerp';
          next unless $project->{species} eq $s_name;
          
          # determine type of data the file has`
          my @types = get_vcf_content_types($project);
          
          # Check if the file have genotype data and being showed
          my $has_sample_data = 1 if ( grep(/^genotype$/, @types) && genotype_samples_exists($s_name, $project, @hostnames) );
          my $has_pop_data = 1 if ( grep(/^populations$/, @types) || is_freq_from_gts($s_name, $project, @hostnames) );
    
          if ( defined $has_sample_data || defined $has_pop_data ){
            # Count the number of variations if the vcf file is used as source
            my $count_var = get_variant_count($project);
            
            # Count difference with previous release
            # TBD - CURRENT PREV DATA FILE HAVE COUNT FOR SOURCE VARIANT DATA 
            # WE NEED SEPARATE GENOTYPE DATA TO HAVE THIS DIFF WORK 
            # my $count_p_var = $prev_data ? $prev_data->{$s_name}->{count_num} : 0;
            $count_p_var = $count_var;
            
            if ($count_var && $count_var > 0){
              # Label
              $species_list{'vcf'}{$s_name}{label} = $label_name;
              
              # Disply name
              $species_list{'vcf'}{$s_name}{'name'} = $display_name;
              $display_list{'vcf'}{$display_name} = $s_name;

              $species_list{'vcf'}{$s_name}{'a'} = $count_var if $has_sample_data;
              $species_list{'vcf'}{$s_name}{'b'} = $count_var if $has_pop_data;
              
              $species_list{'vcf'}{$s_name}{'p_a'} = ($count_var-$count_p_var) if $has_sample_data;
              $species_list{'vcf'}{$s_name}{'p_b'} = ($count_var-$count_p_var) if $has_pop_data;
            }
            
            $b_type = "num";
          }
        }
      }
    }
  }

  # Get html content with derived data and append 
  my ($html_content, $count_species) = generate_html_content(\%species_list,\%display_list,$label_a,$label_b,$b_type);
  $html .= qq{\n  <h2 id="$anchor" style="margin-top:40px">$type data</h2>};
  $html .= q{<p>}.$sql_list{$type}{'extra'}.q{</p>} if ($sql_list{$type}{'extra'});
  $html .= qq{<p style="padding-top:0px;margin-top:0px">There are currently <span style="font-weight:bold;font-size:1.1em;color:#000">$count_species</span> species with $lc_type data in the Ensembl Variation databases:</p>\n};
  $html .= $html_content;  
  
  # Get html content with derived data for vcf species and append 
  if ($sql_list{$type}{'vcf'}){
    my ($html_content, $count_species) = generate_html_content($species_list{'vcf'},$display_list{'vcf'},$label_a,$label_b,$b_type);
    $html .= qq{<p style="padding-top:0px;margin-top:0px">There are currently <span style="font-weight:bold;font-size:1.1em;color:#000">$count_species</span> species with $lc_type data loaded dynamically from vcf file:</p>\n};
    $html .= $html_content;  
  }
}

# Legend
foreach my $type (sort { $colours{$a}{'order'} <=> $colours{$b}{'order'} } keys(%colours)) {
  my $desc  = $colours{$type}{'legend'};
  my $class = $colours{$type}{'colour'};

  $html_legend .= qq{
      <tr>
        <td style="padding-top:4px;text-align:center">
          <span class="vdoc_count_legend $class"></span>
        </td>
        <td style="padding-top:4px">$desc</td>
      </tr>};
}
$html_legend .= qq{    </table>\n  </div>\n</div>};

  

# HTML/output file
open  HTML, "> $html_file" or die "Can't open $html_file : $!";
print HTML $html_header."\n";
print HTML $html_legend."\n";
print HTML $html_title."\n";
print HTML $html."\n";
print HTML qq{<p style="padding-top:15px">The <b>full list of species</b> with their assembly versions in Ensembl is available <a href="/info/about/species.html">here</a>.</p>\n};
print HTML $html_footer."\n";
close(HTML);



###############################################################
##########             FUNCTIONS                     ##########
###############################################################

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
  my $type = shift;
     $type ||= 'variant';
     $type .= 's' if ($type !~ /data$/);
  
  my $count_label;
  my $count_display;
  my $bg_class;
  # From 1 to 9.9 million
  if ($count =~ /^(\d)(\d)\d{5}$/) {
    my $number = ($2!=0) ? "$1.$2" : $1;
    $count = "$number M";
    $count_label = "Over $number million $type";
    $count_display = $count;
    $bg_class = $colours{'few_million'}{'colour'};
  }
  # From 10 tp 99.9 million
  elsif ($count =~ /^(\d{2})\d{6}$/) {
    my $number = $1;
    $count = "$number M";
    $count_label = "Over $number million $type";
    $count_display = $count;
    $bg_class = $colours{'lot_million'}{'colour'};
  }
  # From 100 million
  elsif ($count =~ /^(\d{3}\d*)\d{6}$/) {
    my $number = $1;
    $count = "$number M";
    $count_label = "Over $number million $type";
    $count_display = $count;
    $bg_class = $colours{'hundred_million'}{'colour'};
  }
  # From 1,000 to 999,999
  elsif ($count =~ /^(\d+)\d{3}$/) {
    my $number = $1;
    $count = "$number K";
    $count_label = "Over $number,000 $type";
    $count_display = $count;
    $bg_class = $colours{'thousand'}{'colour'};
  }
  # From 1 to 999
  elsif ($count > 0) {
    $count_label = "$count $type";
    $count_display = "$count";
    $bg_class = $colours{'hundred'}{'colour'};
  }
  # No data
  else {
    $count_label = "No data";
    $count_display = "$count";
    $bg_class = $colours{'zero'}{'colour'};
  }
  return qq{<span class="vdoc_var_count $bg_class" title="$count_label">$count_display</span>};
}

sub round_count_diff {
  my $count = shift;
  my $type  = shift;
     $type ||= 'variant';
     $type .= 's';
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


sub set_bg {
  return ($bg eq '') ? ' class="bg2"' : '';
}


sub generate_html_content {
  my ($species_list, $display_list, $label_a, $label_b, $b_type) = @_;
  
  my $lc_label_a = lc($label_a);
  my $lc_label_b = lc($label_b);
  
  my $html_content = qq{
    <table class="ss" style="width:auto">
      
      <tr class="ss_header">
        <th style="$th_bg">Species</th>  
        <th style="$th_bg;$th_border_left">$label_a</th>
        <th style="$th_bg;padding-left:0px">
          <span class="_ht ht" title="$label_a count difference with the previous Ensembl release (v.$p_version)">
            <small>(e!$p_version &rarr; e!$e_version)</small>
          </span>
        </th>
        <th style="$th_bg;$th_border_left">$label_b</th>
     };
  if ($b_type eq 'num') {
       $html_content .= qq{   
        <th style="$th_bg;padding-left:0px">
          <span class="_ht ht" title="$label_b count difference with the previous Ensembl release (v.$p_version)">
            <small>(e!$p_version &rarr; e!$e_version)</small>
          </span>
        </th>\n};
  }
  $html_content .= qq{      </tr>};
  
  # Loop to display the data by species
  $bg = '';
  my $count_species = 0;
  foreach my $display_name (sort keys(%{ $display_list })) {
    
    next if $display_name =~ /vcf/;

    my $sp = $display_list->{$display_name};
    my $label = $species_list->{$sp}->{'label'};
    my $uc_sp = ucfirst($sp);      
    my $img_src = "/i/species/$uc_sp.png";
    my $img_class = "badge-32";
    my $display_name = $species_list->{$sp}->{'name'};
    my $a_count = $species_list->{$sp}->{'a'};
    my $b_count = $species_list->{$sp}->{'b'};
    
    next if ($a_count eq "0" && ($b_count eq "0" && $b_type eq 'num'));
    next if ($a_count eq "0" && $b_type ne 'num');
    
    $count_species++;
    $a_count = round_count($a_count,$lc_label_a);
    $b_count = ($b_type eq 'num') ? round_count($b_count,$lc_label_b) : qq{<ul style="margin-bottom:0px"><li style="margin-top:0px">}.join("</li><li>", split(',',$b_count))."</li></ul>";
    my $a_p_count = round_count_diff($species_list->{$sp}->{'p_a'},$lc_label_a);
    my $b_p_count = ($b_type eq 'num') ? round_count_diff($species_list->{$sp}->{'p_b'},$lc_label_b) : undef;
  
    my $b_align = ($b_type eq 'num') ? 'right' : 'left';
    if ($b_type ne 'num') {
      $b_count =~ s/StructuralVariation/Structural Variation/;
      $b_count =~ s/SupportingStructuralVariation/Supporting Structural Variation/;
    }
    
    # Species data
    $html_content .= qq{
    <tr$bg style="vertical-align:middle">
      <td>
        <div>
          <div style="float:left">
            <a href="/$uc_sp/Info/Index" title="$label Ensembl Home page" style="vertical-align:middle" target="_blank">
              <img src="$img_src" alt="$label" class="$img_class" style="vertical-align:middle" />
            </a>
           </div>
           <div style="float:left;margin-left:6px;padding-top:2px">
             <div class="bigtext">$display_name</div>
             <div class="small" style="font-style:italic">$label</div>
           </div>
           <div style="clear:both"></div>
         </div>
      </td>
      <td style="text-align:right">$a_count</td>
      <td style="text-align:right">$a_p_count</td>
      <td style="text-align:$b_align">$b_count</td>\n};
    $html_content .= qq{      <td style="text-align:right">$b_p_count</td>\n} if ($b_type eq 'num');
    $html_content .= qq{  </tr>};
    $bg = set_bg();
  }
  $html_content .= qq{</table>\n};

  return ($html_content, $count_species);
}


sub usage {
  
  print qq{
  Usage: perl species_list_detailed_counts.pl [OPTION]
  
  Put detailed data information for each species, into an HTML document.
  
  Options:

    -help           Print this message
      
    -v              Ensembl version, e.g. 65 (Required)
    -o              An HTML output file name (Required)      
    -hlist          The list of host names where the new databases are stored, separated by a coma,
                    e.g. ensembldb.ensembl.org1, ensembldb.ensembl.org2 (Required)
    -phost          Host name where the previous databases are stored, e.g. ensembldb.ensembl.org  (Required)
    -user           MySQL user name (Required)
    -port           MySQL port. 3306 by default (optional)
    -config         The location of the vcf_config.json file. By default it uses the existing one in the same repository.
    -d_dir          The directory location of where the local vcf files are stored (optional). By default it looks in - 
                    /nfs/production/flicek/ensembl/production/ensemblftp/data_files/vertebrates 
    -p_data         Location of the json file that contain result of the previous release for comparison. (Required)
    -site           The URL of the website (optional)
    -etype          The type of Ensembl, e.g. Plant (optional)
  } . "\n";
  exit(0);
}
