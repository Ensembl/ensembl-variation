#!/usr/bin/env perl

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
my ($e_version,$html_file,$hlist,$user,$help,$d_dir);
## EG options
my ($site, $etype);

usage() if (!scalar(@ARGV));
 
GetOptions(
     'v=s'       => \$e_version,
     'o=s'       => \$html_file,
     'help!'     => \$help,
     'hlist=s'   => \$hlist,
     'user=s'    => \$user,
     'd_dir=s'  => \$d_dir
);

## Missing arguments ##
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

# Get the dir this script is residing in
my $dirname = dirname(__FILE__);

# Get the local dir where the vcf files are located
my $data_dir = "/nfs/production/flicek/ensembl/production/ensemblftp/data_files/vertebrates";

if ($d_dir){
  $data_dir = $d_dir;
}

my $vcf_config_file = $dirname . '/../../modules/Bio/EnsEMBL/Variation/DBSQL/vcf_config.json';

# read config from JSON config file
open IN, $vcf_config_file or throw("ERROR: Could not read from config file $vcf_config_file");
local $/ = undef;
my $json_string = <IN>;
 close IN;
    
# parse JSON into hashref $config
my $vcf_config = JSON->new->decode($json_string) or throw("ERROR: Failed to parse config file $vcf_config_file");



## Settings ##
           
my %project_urls = (
  '1000 Genomes'    => 'http://www.1000genomes.org',
  'gnomAD'          => 'http://gnomad.broadinstitute.org/',
  'TOPMed'          => 'https://www.nhlbi.nih.gov/research/resources/nhlbi-precision-medicine-initiative/topmed',
  'UK10K'           => 'https://www.uk10k.org/',
  'MGP'             => 'http://www.sanger.ac.uk/resources/mouse/genomes/',
  'NextGen Project' => 'http://projects.ensembl.org/nextgen/',
  'EVA'             => 'https://www.ebi.ac.uk/eva/?eva-study=###ID###',
  'Gambian Genome Variation Project' => 'https://www.internationalgenome.org/data-portal/data-collection/ggvp-grch38',
  'NCBI ALFA'       => 'https://www.ncbi.nlm.nih.gov/snp/docs/gsr/alfa/',
  'GEM-J'           => 'https://grch38.togovar.org/doc/datasets/gem_j_wga/',
  'NHLBI Exome Sequencing Project' => 'https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000403.v3.p3',
  'ISGC'            => 'http://www.sheephapmap.org/'
);

my $server_name = 'https://static.ensembl.org';
   $server_name = $site if ($site) ;
my $ecaption = 'Ensembl';
my @hostnames = split /,/, $hlist;
my $database = "";
my $pswd = "";
my $db_type = 'variation';

my $margin_bottom_max = '35px';
my $size_max_width = '36px';
my $img_class = "badge-48";

my $evidence_icon_prefix = '/i/val/evidence_';
my $evidence_icon_suffix = '.png';
my $evidence_doc_url  = '../prediction/variant_quality.html#evidence_status';

my $sql  = qq{SHOW DATABASES LIKE '%$db_type\_$e_version%'};
my $sql2 = qq{SELECT p.name,p.description,p.size,p.freqs_from_gts,d.display_name,d.display_priority FROM population p, display_group d WHERE p.display_group_id=d.display_group_id ORDER by d.display_priority};

my $bg = '';
my %pops_list;
my %species_host;
my %species_usual_name;
my %species_subpop;
my @evidence_list;


## Headers ##
my $pop_table_header = qq{
  <tr>
    <th style="width:200px">Name</th>
    <th style="width:$size_max_width">Size</th>
    <th>Description</th>
  </tr>
};


###############################################################
##########             MAIN PART                       ########
###############################################################

## Species / host / database ##
foreach my $hostname (@hostnames) {
  my $sth = get_connection_and_query($database, $hostname, $sql);

  # loop over databases
  while (my ($dbname) = $sth->fetchrow_array) {
    next if ($dbname !~ /^[a-z][a-z_]*_[a-z]+_$db_type\_$e_version\_\d+$/i);
    next if ($dbname =~ /^(master_schema|drosophila|saccharomyces|ciona)/ || $dbname =~ /^homo_sapiens_$db_type\_\d+_37$/ || $dbname =~ /private/);

    print STDERR "$dbname\n";
    $dbname =~ /^(.+)_$db_type/;
    my $s_name = $1; 
    
    my $label_name = ucfirst($s_name);
       $label_name =~ s/_/ /g;
    
    # Get list of triplets: species / DB / Host
    $species_host{$label_name} = {'host' => $hostname, 'dbname' => $dbname};

    # Get the sub populations
    $species_subpop{$label_name} = get_sub_populations($label_name);

    # Get the list of evidence status
    if (!@evidence_list) {
      @evidence_list = @{get_evidence($dbname,$hostname)};
    }
  }
  $sth->finish;
}

## Get all the data in (organised) hashes
get_project_populations();


## Populations ##
my $html_pop = '';


# Loop over the species (placing human first)
foreach my $species (sort { ($a !~ /Homo/ cmp $b !~ /Homo/) || $a cmp $b } keys(%pops_list)) {
  next unless  %{ $pops_list{$species} };

  my $id_species = $species;
     $id_species =~ s/ /_/g;
  my $species_label = $species_usual_name{$species};

  my $margin_top = ($species =~ /Homo/i) ? '20px' : '50px';

  # Species header
  $html_pop .= qq{
    <div style="padding-left:0px;padding-bottom:3px;margin-top:$margin_top;margin-bottom:20px">
      <a href="/$id_species/Info/Index" title="$species Ensembl Home page" style="vertical-align:middle" target="_blank"><img src="/i/species/$id_species.png" alt="$id_species" class="$img_class" style="float:none;margin-right:4px;vertical-align:middle" /></a>
      <h2 id="$id_species" style="display:inline;color:#333">$species<span class="small vdoc_species_sci_name"> ($species_label)</span></h2>
    </div>
    <div style="margin-left:10px">
  };

  # Loop over the projects (from the VCF config file)
  foreach my $project_label (sort { $a cmp $b } keys(%{$pops_list{$species}})) {

    my %pop_seen;

    my $project_id = lc($project_label);
    $project_id =~ s/ /_/g;
    my $html_current_pop = qq{<table id="$project_id" class="ss" style="margin-bottom:5px">\n  $pop_table_header\n};

    my $evidence = $pops_list{$species}{$project_label}{'evidence'};

    my $count_project_entries = 0;

    # Loop over the different files/parts of a same project
    foreach my $project (sort { $a cmp $b } keys(%{$pops_list{$species}{$project_label}})) {

      $bg = '';

      next if (!$pops_list{$species}{$project_label}{$project}{'pop_list'});

      my %pop_tree;
      my %sub_pops;

      # Get the populations structure (if it exists)
      %pop_tree = %{$pops_list{$species}{$project_label}{$project}{'pop_tree'}} if ($pops_list{$species}{$project_label}{$project}{'pop_tree'});
      %sub_pops = %{$pops_list{$species}{$project_label}{$project}{'sub_pops'}} if ($pops_list{$species}{$project_label}{$project}{'sub_pops'});

      my %pop_data = %{$pops_list{$species}{$project_label}{$project}{'pop_data'}};
      my @pop_list = @{$pops_list{$species}{$project_label}{$project}{'pop_list'}};

      @pop_list = sort { ($a !~ /ALL/ cmp $b !~ /ALL/) || $a cmp $b } @pop_list if (!%pop_tree);

      # Loop over the populations and add a row for each
      foreach my $pop (@pop_list) {
        next if ($pop_seen{$pop});

        # Avoid duplicated entries within the same project
        $pop_seen{$pop} = 1;

        my $new_bg = $bg;
        if ($pop_tree{$pop}) {
          $new_bg = ($pop =~ /all$/i) ? ' class="supergroup"' : ' class="subgroup"';
        }
        my $p_name = $pop_data{$pop}{'label'};
        if (!$pop_tree{$pop} && %pop_tree && $sub_pops{$pop}) {
           $p_name = qq{<ul style="margin:0px"><li style="margin:0px">$p_name</li></ul>};
        }
        my $desc   = $pop_data{$pop}{'desc'};
        my $size   = $pop_data{$pop}{'size'};

        # Create a HTML row with the population data
        $html_current_pop .= qq{  <tr$new_bg>\n    <td>$p_name</td>\n    <td style="text-align:right">$size</td>\n    <td>$desc</td>\n  </tr>\n};
        $bg = set_bg($bg);
        $count_project_entries ++;
      }
    }
    $html_current_pop .= "</table>\n";

    # Create the project header
    my $plural = ($count_project_entries > 1) ? 's' : '';

    my $url;
    foreach my $project_name (keys(%project_urls)) {
      if ($project_label =~ /$project_name/) {
        $url = $project_urls{$project_name};
        if ($url =~ /###ID###/ && $project_label =~ /EVA study\s(\w+)$/) {
          my $study_id = $1;
          $url =~ s/###ID###/$study_id/;
        }
        last;
      }
    }
    my $project_title = ($url) ? qq{<a href="$url" target="_blank" style="text-decoration:none">$project_label</a>} : $project_label;
       $project_title = ' the '.$project_title if ($project_title =~ /(project|study|consortium)/i);

    $html_pop .= qq{
      <h3>Population$plural from $project_title</h3>
      $html_current_pop\n};


    # Evidence status
    my $has_evidence = 0;
    foreach my $evidence (@evidence_list) {
      my $label_no_space = $project_label;
         $label_no_space =~ s/ //g;

      if ($project_label =~ /$evidence/ || $label_no_space =~ /$evidence/) {
        my $evidence_img = "$evidence_icon_prefix$evidence$evidence_icon_suffix";
        $html_pop .= qq{
      <p style="margin-bottom:$margin_bottom_max">
      Variants which have been discovered in this project have the "evidence status" <a href="$evidence_doc_url"><b>$evidence</b></a>.
      On the website this corresponds to the icon <a href="$evidence_doc_url"><img class="_ht" src="$evidence_img" title="$evidence" style="vertical-align:bottom"/></a>.
      </p>};
      $has_evidence = 1;
      last;
      }
    }

    if ($has_evidence == 0) {
      $html_pop .= qq{\n  <div style="margin-bottom:$margin_bottom_max"></div>\n};
    }
  }
  $html_pop .= qq{</div>};
}


## HTML/output file ##
open  HTML, "> $html_file" or die "Can't open $html_file : $!";
print HTML $html_pop;
close(HTML);




###############################################################
##########             FUNCTIONS                     ##########
###############################################################


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

# Determine what type data contains in the vcf file
sub get_vcf_content_types {
  my ($project) = @_;
  my @types;
  
  # this ignores the false positive sigpipe error from tabix command 
  $SIG{PIPE} = 'DEFAULT';
  
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

  my $genotypes = `tabix -D $file_full_path -H | grep '##FORMAT' | grep 'ID=GT'`;
  push @types, "genotype" if $genotypes;
  
  # check in a actual line for FORMAT field if not exist in header
  unless ($genotypes){
    my $chr = `tabix -D $file_full_path -l | head -n 1`;
    chop $chr;
  
    my $line = `tabix -D $file_full_path $chr | head -n 1`;
  
    my $format_field = (split /\t/, $line)[8];
    push @types, "genotype" if $format_field;

    my $info_field = (split /\t/, $line)[7];
    if ( ($info_field =~ /AF=/) || ($info_field =~ /AC=/ && $info_field =~ /AN=/) ) {
      push @types, "frequency";
    }
    # a hard-coded check for NCBI-ALPHA and TOPMED as they have very special field for frequency
    if ( ($info_field =~ /AN_SAMN/) || ($info_field =~ /TOPMED=/) ) {
      push @types, "frequency";
    }
  }
  
  return @types;
}

# Build the project populations structure if it exists
sub get_population_structure {
  my $pops     = shift;
  my $pop_tree = shift;
  my $sub_pops = shift;

  my $pop_list;
  if ($pop_tree) {

    foreach my $pop_id (sort { ($a !~ /ALL/ cmp $b !~ /ALL/) || $a cmp $b } keys %$pops) {
      next if (grep { $_ eq $pop_id} @$pop_list);

      if (!$pop_tree->{$pop_id} && !$sub_pops->{$pop_id}) {
        push (@$pop_list, $pop_id);
      }
      elsif ($pop_tree->{$pop_id}) {
        push (@$pop_list, $pop_id);
        $pop_list = add_population_to_list($pop_tree,$pop_id,$pop_list);
      }
    }
  }
  else {
    my @list = sort { ($a =~ /ALL/ cmp $b =~ /ALL/) || $a cmp $b } keys %$pops;
    $pop_list = \@list;
  }
  return $pop_list;
}

# Add the sub-populations after the super-population in the "pop_list" array
sub add_population_to_list {
  my $pop_tree = shift;
  my $pop_id   = shift;
  my $pop_list = shift;

  foreach my $sub_pop_id (sort {$a cmp $b} keys(%{$pop_tree->{$pop_id}})) {
    push (@$pop_list, $sub_pop_id);
    if ($pop_tree->{$sub_pop_id}) {
      $pop_list = add_population_to_list($pop_tree,$sub_pop_id,$pop_list);
    }
  }
  return $pop_list;
}

sub get_evidence {
  my $dbname   = shift;
  my $hostname = shift;
  
  my @evidence_list;
  
  my $stmt = qq{ SELECT value FROM attrib WHERE attrib_type_id=497 };
  my $sth = get_connection_and_query($dbname, $hostname, $stmt);

  while(my @data = $sth->fetchrow_array) {
    push @evidence_list, $data[0];
  }
  $sth->finish;
  return \@evidence_list;
}


# Get population size if it hasn't been set in the "population" SQL table
sub get_size {
  my $pop_id   = shift;
  my $dbname   = shift;
  my $hostname = shift;
 
  my $stmt = qq{ SELECT count(*) FROM sample_population WHERE population_id=?};
  my $sth = get_connection_and_query($dbname, $hostname, $stmt, [$pop_id]);
  my $size = ($sth->fetchrow_array)[0];
  $sth->finish;

  return ($size == 0) ? '-' : $size;
}

# Alternate the background colour of the table rows
sub set_bg {
  my $bg_class = shift;
  return ($bg_class eq '') ? ' class="bg2"' : '';
}


# Format the population description text
sub parse_desc {
  my $content = shift;
  $content = '-' if (!$content);
  my @desc = split(/\.,/, $content);
  $content = "$desc[0]. $desc[1]." if scalar(@desc > 1);
  
  return $content;
}

# Store the sub population names 
sub get_sub_populations {
  my $species = shift;
  
  my $host   = $species_host{$species}{'host'};
  my $dbname = $species_host{$species}{'dbname'};
  
  my %subpop_data;
  
  my $pop_sql = qq{ SELECT p1.population_id, p2.name FROM population p1, population p2, population_structure ps 
                    WHERE p1.population_id=ps.super_population_id AND p2.population_id=ps.sub_population_id };
  my $pop_sth = get_connection_and_query($dbname, $host, $pop_sql);
  while(my @data = $pop_sth->fetchrow_array) {
    $subpop_data{$data[0]}{$data[1]} = 1;
  }
  $pop_sth->finish;
  return \%subpop_data;
}

# Build a hash containing all the relevant information to create the species/projects/populations content
# Mostly based on the vcf_config file and Variation databases
sub get_project_populations {

  foreach my $project (@{$vcf_config->{'collections'}}) {
    # Check if the file have genotype data and being showed
    my @types = get_vcf_content_types($project);
    next unless ( grep(/^genotype$/, @types) || grep(/^frequency$/, @types) );
    
    my $project_id = $project->{'id'};
    next if ($project->{'assembly'} =~ /GRCh37/i || $project->{'annotation_type'} eq 'cadd' || $project->{'annotation_type'} eq 'gerp');

    my $species = ucfirst($project->{'species'});
       $species =~ s/_/ /g;

    if (!$species_usual_name{$species}) {
      my $spe_host   = $species_host{$species}{'host'};
      my $spe_dbname = $species_host{$species}{'dbname'};
         $spe_dbname =~ s/variation/core/;
      my $spe_stmt = qq{ SELECT meta_value FROM meta WHERE meta_key='species.display_name'};
 
      my $spe_sth  = get_connection_and_query($spe_dbname, $spe_host, $spe_stmt);
      $species_usual_name{$species} = ($spe_sth->fetchrow_array)[0];
      $spe_sth->finish;
    }

    my $project_label = get_project_label($project,$species);

    my $population_prefix;
    if ($project->{'population_prefix'}) {
      $population_prefix = $project->{'population_prefix'};
    }
    elsif ($project->{'sample_prefix'}) {
      $population_prefix = $project->{'sample_prefix'};
    }
    # Special case for MGP
    if ($population_prefix && $population_prefix =~ /^MGP:/) {
      $population_prefix = "Mouse Genomes Project";
    }

    my %pop_data;
    my %pop_tree;
    my %sub_pops;
    my $pop_list;

    if ($project->{'populations'}) {
      foreach my $pop (keys(%{$project->{'populations'}})) {
        next if ($pop !~ /\w+/ || $pop eq '');
        my $pop_name = $project->{'populations'}{$pop}{'name'};
        if ($population_prefix && $pop_name !~ /$population_prefix/i) {
          $pop_name = $population_prefix.$pop_name;
        }
        my @composed_name = split(':', $pop_name);
           $composed_name[$#composed_name] = '<b>'.$composed_name[$#composed_name].'</b>';
        $pop_name = join(':',@composed_name);

        my $pop_desc = $project->{'populations'}{$pop}{'description'};

        $pop_data{$pop_name} = {'id'    => $pop,
                                'label' => $pop_name,
                                'desc'  => $pop_desc,
                                'size'  => '-'
                               };
        push(@$pop_list,$pop_name);
      }
    }
    else {
      my $term = ($population_prefix) ? $population_prefix : '';
         $term =~ s/:$//;
      my $dbname = $species_host{$species}{'dbname'};
      my $host   = $species_host{$species}{'host'};
      my $stmt = qq{ SELECT population_id, name, size, description, display_group_id FROM population WHERE name like ? or name = ? ORDER BY name};

      my $sth  = get_connection_and_query($dbname, $host, $stmt, ["$population_prefix%",$term]);

      while(my @data = $sth->fetchrow_array) {
        # Skip if the population is not in display group
        next unless defined $data[4];

        my @composed_name = split(':', $data[1]);
           $composed_name[$#composed_name] = '<b>'.$composed_name[$#composed_name].'</b>';
        my $pop_name = join(':',@composed_name);

        $data[2] = '-' if (!$data[2]);
        my $desc = parse_desc($data[3]);
        my $size = ($data[2] && $data[2] ne '-' ) ? $data[2] : get_size($data[0], $dbname, $host);

        $pop_data{$data[1]} = {'id'    => $data[0],
                               'label' => $pop_name,
                               'desc'  => $desc,
                               'size'  => $size
                              };
        # Super/sub populations   
        if ($species_subpop{$species}{$data[0]}) {
          foreach my $sub_pop (keys(%{$species_subpop{$species}{$data[0]}})) {
            if (($population_prefix && $sub_pop =~ /^$population_prefix/) || !$population_prefix) {
              $sub_pops{$sub_pop} = 1;
              $pop_tree{$data[1]}{$sub_pop} = 1;
            }
          }
        }
      }
      $sth->finish;
      $pop_list = get_population_structure(\%pop_data, \%pop_tree, \%sub_pops);
    }

    $pops_list{$species}{$project_label}{$project}{'pop_data'} = \%pop_data if %pop_data ;
    $pops_list{$species}{$project_label}{$project}{'pop_tree'} = \%pop_tree if %pop_tree;
    $pops_list{$species}{$project_label}{$project}{'sub_pops'} = \%sub_pops if %sub_pops;
    $pops_list{$species}{$project_label}{$project}{'pop_list'} = $pop_list if defined $pop_list;
  }
}

# Get a more displayable project label from the project ID/name in the vcf_config file
sub get_project_label {
  my $project = shift;
  my $species = shift;

  my $label =  $project->{'id'};
     $label =~ s/_GRCh38//i;

  if ($project->{'population_display_group'} && $project->{'population_display_group'}{'display_group_name'}) {
    $label = $project->{'population_display_group'}{'display_group_name'};
  }

  if ($label =~ /^1000/) {
    $label = '1000 Genomes Project';
  }
  elsif ($label =~ /^nextgen/) {
    $label = 'NextGen Project';
  }
  elsif ($label =~ /EVA_(.+)$/) {
    $label = "EVA study $1";
  }
  elsif ($label =~ /^PRJEB(.+)/) {
    $label = "EVA study $label";
  }
  elsif ($label =~ /mouse_genome_project/) {
    $label = "Mouse Genomes Project (MGP)";
  }
  elsif ($label =~ /sheep_genome_consortium/) {
    $label = "International Sheep Genome Consortium (ISGC)";
  }
  elsif ($label =~ /gambian_genome_variation_project/) {
    $label = "Gambian Genome Variation Project";
  }
  elsif ($project->{'source_name'}) {
    $label = $project->{'source_name'};
  }
  return $label;
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
  my $dbh = DBI->connect($dsn, $user, $pswd) or die "Connection failed";

  my $sth = $dbh->prepare($sql);
  if ($params) {
    $sth->execute(@$params);
  }
  else {
    $sth->execute;
  }
  return $sth;
}


sub usage {

  print qq{
  Usage: perl sources2html.pl [OPTION]

  Create HTML tables for listing the population in the main genotyping projects available in Ensembl Variation.

  Options:

    -help           Print this message
    -v              Ensembl version, e.g. 65 (Required)
    -o              An HTML output file name (Required)      
    -hlist          The list of host names (with port) where the new databases are stored, separated by a coma,
                    e.g. ensembldb.ensembl.org1:1234, ensembldb.ensembl.org2:1234 (Required)
    -user           MySQL user name (Required)
    -d_dir          The directory location of where the local vcf files are stored (optional). By default it looks in - 
                    /nfs/production/flicek/ensembl/production/ensemblftp/data_files/vertebrates 
  } . "\n";
  exit(0);
}

