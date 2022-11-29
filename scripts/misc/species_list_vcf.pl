# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2022] EMBL-European Bioinformatics Institute
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

###############################################################
##########             CONFIGURE                        #######
###############################################################

my ($e_version, $html_file, $config, $d_dir, $p_data, $hlist, $user, $dump_file, $help);
## EG options
my ($site, $etype);

usage() if (!scalar(@ARGV));
 
GetOptions(
     'v=s'      => \$e_version,
     'o=s'      => \$html_file,
     'dump=s'   => \$dump_file,
     'config=s' => \$config,
     'd_dir=s'  => \$d_dir,
     'p_data=s' => \$p_data,
     'hlist=s'  => \$hlist,
     'user=s'   => \$user,
     'site=s'   => \$site,
     'etype=s'  => \$etype,
     'help!'    => \$help
);

if (!$e_version) {
  print "> Error! Please give an Ensembl version, using the option '-v' \n";
  usage();
}
if (!$html_file) {
  print "> Error! Please give an output file using the option '-o'\n";
  usage();
}
if (!$dump_file) {
  print "> Error! Please give an dump file using the option '-dump'\n";
  usage();
}
if (!$p_data) {
  print "> Error! Please give the location of previous release data using the option '-p_data'\n";
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

if ($site) {
  $server_name = $site;
}

my @hostnames = split /,/, $hlist;

# Get the vcf config file location
my $dirname = dirname(__FILE__);
my $vcf_config_file = $dirname . '/../../modules/Bio/EnsEMBL/Variation/DBSQL/vcf_config.json';

if ($config){
  $vcf_config_file = $config;
}

# Get the local dir where the vcf files are located
my $data_dir = "/nfs/production/flicek/ensembl/production/ensemblftp/data_files/vertebrates";

if ($d_dir){
  $data_dir = $d_dir;
}

my $pswd = "";
my $p_version = $e_version-1;
my $detailed_counts = 'species_detailed_counts.html';

my $html;
   
my %colours = ( 'hundred_million' => { 'order' => 5, 'colour' => 'vdoc_million_0', 'legend' => 'From 100 million'},
                'lot_million'     => { 'order' => 4, 'colour' => 'vdoc_million_1', 'legend' => 'From 10 million to 99,9 million'},
                'few_million'     => { 'order' => 3, 'colour' => 'vdoc_million_2', 'legend' => 'From 1 million to 9,9 million'},
                'thousand'        => { 'order' => 2, 'colour' => 'vdoc_thousand',  'legend' => 'From 1,000 to 999,999'},
                'hundred'         => { 'order' => 1, 'colour' => 'vdoc_hundred',   'legend' => 'From 1 to 999'}
              );              
              
my %tables = ( 'Sample'             => { 'order' => 2 , 'anchor' => '#genotype',             'table' => 'compressed_genotype_var'},
               'Population'         => { 'order' => 3 , 'anchor' => '#genotype',             'table' => 'population_genotype'},
               'Phenotype'          => { 'order' => 4 , 'anchor' => '#phenotype',            'table' => 'phenotype_feature'},
               'Citation'           => { 'order' => 5 , 'anchor' => '#citation',             'table' => 'variation_citation'},
               'Structural variant' => { 'order' => 1 , 'anchor' => '#structural_variation', 'table' => 'structural_variation'}
             );
                        
my %species_list;
my %display_list;

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

  my $prev_data = JSON->new->decode($json_string) or throw("ERROR: Failed to parse config file $p_data");
}

###############################################################
##########             MAIN PART                       ########
###############################################################

foreach my $project (@{ $vcf_config->{'collections'} }) {
  next if $project->{annotation_type} eq 'cadd' || $project->{annotation_type} eq 'gerp';

  my $s_name = $project->{species};

  if ($etype) { # EG site - need to filter out species
    my $img_thumb = sprintf qq{eg-plugins/%s/htdocs/img/species/thumb_%s.png}, $etype, ucfirst($s_name);
    #  print "- checking for $img_thumb ... ";
    if (! -e $img_thumb) {
      print "\t... skipping \n";
      next;
    } 
  }

  # determine type of data the file has
  my @types = get_vcf_content_types($project);
    
  # We are only interested with species which are vcf-only for now
  if ( grep /^source$/, @types){
    # Count the number of variations if the vcf file is used as source
    my $count_var = get_variant_count($project);
    if ($count_var && $count_var > 0){
      $species_list{$s_name}{count} = round_count($count_var);
    }

    # Check if the file have genotype data and being showed
    if ( grep /^genotype$/, @types){
      # Check if either vcf config or database have the samples
      if ( genotype_samples_exists($s_name, $project) ){
        $species_list{$s_name}{genotype} = 1;
      }
    }

    # Get the species labels
    my $label_name = ucfirst($s_name);
        $label_name =~ s/_gca[0-9]{9}[v0-9]*+$//g;	# remove any gca from name 
        $label_name =~ s/_/ /g;
    $species_list{$s_name}{label} = $label_name;
  
    # Get species display name
    foreach my $hostname (@hostnames) {
      my $display_name;

      my $sql = qq{SHOW DATABASES LIKE '%$s_name\%core\_$e_version%'};
      my $sth = get_connection_and_query("", $hostname, $sql);

      while (my ($core_dbname) = $sth->fetchrow_array) {
        my $sql2 = qq{SELECT meta_value FROM meta WHERE meta_key="species.display_name" LIMIT 1};
        my $sth2 = get_connection_and_query($core_dbname, $hostname, $sql2);
        $display_name = $sth2->fetchrow_array;  
        
        $species_list{$s_name}{'name'} = $display_name;
        $display_list{$display_name} = $s_name;

        last if $display_name;
      }

      last if $display_name;
    }

    # Count difference with previous release
    if ($prev_data) {
      # count the difference
      my $count_p_var = $prev_data->{$s_name}->{count};
      $count_p_var =~ s/<[^>]*.//g;
      $species_list{$s_name}{'p_count'} = round_count_diff($count_var-$count_p_var);
    }
  }
}

my $count_species = scalar(keys(%species_list));

my $th_bg = qq{background-color:#BBB};
my $th_border_left = qq{border-left:1px solid #DDD};
my $th_border_left_top = qq{style="$th_border_left;text-align:center"};

my $html_content = qq{
  <table class="ss" style="width:auto">
    <tr class="ss_header">
      <th style="$th_bg">Species</th>  
      <th style="$th_bg">Variant</th>
      <th style="$th_bg;padding-left:0px">
        <span class="_ht ht" title="Variant count difference with the previous Ensembl release (v.$p_version)">
          <small>(e!$p_version &rarr; e!$e_version)</small>
        </span>
      </th>
      <th style="$th_bg">
        <a class="_ht" style="text-decoration:none" title="See detailed counts" href=species_detailed_counts.html#genotype>
          Genotype
        </a>
      </th>
    </tr>
};
my $bg = '';

foreach my $display_name (sort keys(%display_list)) {

  my $sp = $display_list{$display_name};
  my $label = $species_list{$sp}{'label'};
  my $uc_sp = ucfirst($sp);      
  my $img_src = "/i/species/$uc_sp.png";
  my $img_class = "badge-32";
  my $display_name = $species_list{$sp}{'name'};
  my $var_count = $species_list{$sp}{'count'};
  my $var_p_count = $species_list{$sp}{'p_count'};

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
    <td style="text-align:right">$var_count</td>
    <td style="text-align:right">$var_p_count</td>\n};
  
  
  my $has_data = ($species_list{$sp}{genotype}) ? qq{<img src="/i/16/check.png" title="Data available" />} : '-';
  $html_content .= qq{    <td style="text-align:center">$has_data</td>\n};
  
  $html_content .= qq{  </tr>};
  $bg = set_bg();
}
$html_content .= qq{</table>\n};


# Legend
my $html_legend = qq{
<div>
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
</div>
};


## HTML/output file ##
open  HTML, "> $html_file" or die "Can't open $html_file : $!\n";
print HTML qq{<p style="padding-top:0px;margin-top:0px">We currently have <span style="font-weight:bold;font-size:1.1em;color:#000">$count_species</span> species supported this way:</p>\n};
print HTML $html_content;
print HTML $html_legend;
close HTML;

# Dump the data object - to be stored in repo and used in next release
my $data_dump = JSON->new->encode(\%species_list) or throw("ERROR: Failed to encode data file for dumping\n");
open DATA_DUMP, "> $dump_file";
print DATA_DUMP $data_dump;
close DATA_DUMP;


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

# Determine what type data contains in the vcf file
sub get_vcf_content_types {
  my ($project) = @_;
  my @types;

  # add if the vcf collection mentions annotation type
  push @types, $project->{annotation_type} if $project->{annotation_type};

  # if use_as_source is set then it is the main source for tracks
  push @types, "source" if $project->{use_as_source};

  # check FORMAT field of the vcf file to see if it has genotype
  my $file = get_random_file($project);

  my $file_full_path = $file;
  if ($project->{type} eq "local"){
    $file_full_path = $data_dir . $file_full_path;
  }

  my $genotypes = `tabix $file_full_path -H | grep '##FORMAT' | grep 'ID=GT'`;
  push @types, "genotype" if $genotypes;
  
  # check in a actual line for FORMAT field if not exist in header
  unless ($genotypes){
    my $chr = `tabix $file_full_path -l | head -n 1`;
    chop $chr;
    
    my $line = `tabix $file_full_path $chr | head -n 1`;

    my $format_field = (split /\t/, $line)[8];
    
    push @types, "genotype" if $format_field;
  }
  
  return @types;
}

# Check if samples from a vcf files exist in either vcf config or database
sub genotype_samples_exists {
  my ($species, $project) = @_;
  
  # Get samples from vcf file
  my @samples;
  foreach my $file (get_all_files($project)){
    push @samples, (split / /, `bcftools query -l $file | xargs | tr -d '\n'`);
  }

  # Get samples from vcf config 
  my @samples_in_vcf = keys %{ $project->{sample_populations} };

  # Check if any of the sample matches
  foreach ( @samples_in_vcf ){
    if (grep /^$_$/, @samples){
      return 1;
    }
  }

  my $samples_str = join ",", (map { "'$_'" } @samples);
  foreach my $hostname (@hostnames) {
    my $sql = qq{SHOW DATABASES LIKE '%$species\%variation\_$e_version%'};
    my $sth = get_connection_and_query("", $hostname, $sql);

    while (my ($var_dbname) = $sth->fetchrow_array) {
      my $sql2 = qq{SELECT name FROM sample WHERE name IN ($samples_str) LIMIT 1};
      my $sth2 = get_connection_and_query($var_dbname, $hostname, $sql2);
      
      return 1 if $sth2->fetchrow_array;
    }
  }

  return 0;
}

# Get number of variant from a vcf file
sub get_variant_count {
  my ($project) = @_;
  my $count;

  foreach my $file (get_all_files($project)){
    $count += `bcftools stats $file | grep -ve '^#' | grep -e 'number of records' | cut -d\$'\t' -f 4`;

    unless($count) {
      # If file is remote try downloading it and count the line number
      if ($file =~ /^http/ || $file =~ /^ftp/) {
        `wget -O vcf.gz $file`;
        $count = `bgzip -d -c vcf.gz | grep -v '^#' | wc -l`;
        `rm vcf.gz`;
      }
      # If file is local no need for downloadin; just count line number 
      else{
        $count = `bgzip -d -c $file | grep -v '^#' | wc -l`;
      }
    }
  }

  return $count;
}

# Get a random file from filename template in vcf collection
sub get_random_file {
  my ($project) = @_;
  my $file;

  my $filename_template = $project->{filename_template};

  if ($filename_template =~ /###CHR###/){
    my $chromosomes = $project->{chromosomes};

    return undef unless $chromosomes;

    my $chr = @{ $chromosomes }[0];
    
    $file = $filename_template =~ s/###CHR###/$chr/gr;
  }
  else{
    $file = $filename_template
  }

  return $file;
}

# Get all files from filename template in vcf collection
sub get_all_files {
  my ($project) = @_;
  my @files;

  my $filename_template = $project->{filename_template};

  if ($filename_template =~ /###CHR###/){
    my $chromosomes = $project->{chromosomes};

    return undef unless $chromosomes;

    foreach my $chr (@{ $chromosomes }){
      my $file = $filename_template =~ s/###CHR###/$chr/gr;
      push @files, $file;
    }
  }
  else{
    push @files, $filename_template;
  }

  return @files;
}

sub round_count {
  my $count = shift;
  my $type = 'variants';
  
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
  # From 10 to 99.9 million
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
    -dump           An output data dump file name (json format) (Required)  
    -config         The location of the vcf_config.json file. By default it uses the existing one in the same repository.
    -d_dir          The directory location of where the local vcf files are stored (optional). By default it looks in - 
                    /nfs/production/flicek/ensembl/production/ensemblftp/data_files/vertebrates 
    -p_data         Location of the json file that contain result of the previous release for comparison. (Required)
    -hlist          The list of host names (with port) where the new databases are stored, separated by a coma,
                    e.g. ensembldb.ensembl.org1:1234, ensembldb.ensembl.org2:1234 (Required)
    -user           MySQL user name (Required)
    -site           The URL of the website (optional)
    -etype          The type of Ensembl, e.g. Plant (optional)
  } . "\n";
  exit(0);
}