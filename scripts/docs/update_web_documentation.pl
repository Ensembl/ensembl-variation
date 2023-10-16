#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2023] EMBL-European Bioinformatics Institute
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
  <http://www.ensembl.org/Help/Contact>.

=cut


# Script to update the documentation page "Data description".

use strict;
use warnings;
use DBI;
use Getopt::Long;
use File::Basename;

# Print the usage instructions if run without parameters
usage() unless (scalar(@ARGV));

my ($version,$input_dir,$output_dir,$help,$host,$hlist,$user,$port,$phost,$ohost,$no_subdir,$species,$dump_file,$config,$d_dir,$p_data,$sections);

GetOptions(
  'v=i'         => \$version,
  'i=s'         => \$input_dir,
  'o=s'         => \$output_dir,
  'help!'       => \$help,
  'host=s'      => \$host,
  'hlist=s'     => \$hlist,
  'user=s'      => \$user,
  'port=i'      => \$port,
  'phost=s'     => \$phost,
  'ohost=s'     => \$ohost,
  'no_subdir!'  => \$no_subdir,
  'species|s=s' => \$species,
  'dump=s'      => \$dump_file,
  'config=s'    => \$config,
  'd_dir=s'     => \$d_dir,
  'p_data=s'    => \$p_data,
  'sections=s'  => \$sections
);

usage("input and output directories must be specified") unless ($input_dir && $output_dir);
usage("Host, port and version must be specified") unless ($host && $port && $version);
usage("Hosts list, user must be specified") unless ($hlist && $user);
usage("Previous host must be specified") unless ($phost);
usage("Host providing the ontology database must be specified") unless ($ohost);

# Check if output directory exists
die "Could not find output dir $output_dir - please create it first\n" unless -d $output_dir;

# Get the dir that this script is residing in
my $dirname = dirname(__FILE__);

# Get the vcf config file 
my $vcf_config_file = $dirname . '/../../modules/Bio/EnsEMBL/Variation/DBSQL/vcf_config.json';
if ($config){
  $vcf_config_file = $config;
}

# Get the local dir where the vcf files are located
my $data_dir = "/nfs/production/flicek/ensembl/production/ensemblftp/data_files/vertebrates";
if ($d_dir){
  $data_dir = $d_dir;
}

$species ||= 'Homo_sapiens';
$dump_file ||= "vcf_species_list_${version}.json";

my $prev_version = $version - 1;
$p_data ||= $dirname . "/data_files/vcf_species_list_${prev_version}.json";

$sections ||= "all";
my @sections = split /,/, $sections;

my ($section,$tmp_file,$tmp_section, $file_name, $subdir,$copy2subdir,$content_before,$new_content,$content_after);
my @ontologies = ('cmo','efo','hp','mp','vt');
my %subdirs = ( 'species_data_types.html'              => 'species',
                'sets.html'                            => 'species',
                'populations.html'                     => 'species',
                'sources_documentation.html'           => 'species',
                'classification.html'                  => 'prediction',
                'phenotype_annotation.html'            => 'phenotype',
                'sources_phenotype_documentation.html' => 'phenotype',
                'species_detailed_counts.html'         => 'species'
              );

#### Generates the "List of database species" table documentation
$section = 'species';

if ( $sections eq "all" || grep(/^${section}$/, @sections) ){
  $tmp_file    = "data_desc_$section.html";
  $tmp_section = "$section\_tmp.html";
  $file_name   = "species_data_types.html";
  $subdir      = $subdirs{$file_name};

  print localtime() . "\t# Start species database list...\n";
  `cp $input_dir/$subdir/$file_name $tmp_file`;
  $content_before = get_content($section,'start');
  $content_after  = get_content($section,'end');
  `perl $dirname/species_list.pl -v $version -o $tmp_section -hlist $hlist -user $user -phost $phost`;
  $new_content = `cat $tmp_section`;
  `rm -f $tmp_section`;
  print_into_tmp_file($tmp_file,$content_before,$new_content,$content_after);

  $copy2subdir = ($no_subdir) ? '' : $subdirs{$file_name};
  copy_updated_file($copy2subdir,$file_name,$tmp_file);

  print localtime() . "\t\t> Species database list - finished\n";
}


#### Generates the "Dynamically loaded species" table documentation
$section = 'vcf_species';

if ( $sections eq "all" || grep(/^${section}$/, @sections) ){
  $tmp_file    = "data_desc_$section.html";
  $tmp_section = "$section\_tmp.html";
  $file_name   = "species_data_types.html";
  $subdir      = $subdirs{$file_name};

  print localtime() . "\t# Start species vcf-only list...\n";
  # Use output file from the list of database species update
  if (-d "$output_dir/$subdir") {
    `cp $output_dir/$subdir/$file_name $tmp_file`;
  }
  elsif (-f "$output_dir/$file_name") {
    `cp $output_dir/$file_name $tmp_file`;
  }
  else {
    `cp $input_dir/$subdir/$file_name $tmp_file`;
  }
  $content_before = get_content($section,'start');
  $content_after  = get_content($section,'end');
  `perl $dirname/species_list_vcf.pl -v $version -o $tmp_section -dump $dump_file -p_data $p_data -config $vcf_config_file -d_dir $data_dir -hlist $hlist -user $user`;
  $new_content = `cat $tmp_section`;
  `rm -f $tmp_section`;
  print_into_tmp_file($tmp_file,$content_before,$new_content,$content_after);

  $copy2subdir = ($no_subdir) ? '' : $subdirs{$file_name};
  copy_updated_file($copy2subdir,$file_name,$tmp_file);

  print localtime() . "\t\t> Species vcf-only list - finished\n";
}

#### Generates the "Variation classes" table documentation
$section = 'classes';

if ( $sections eq "all" || grep(/^${section}$/, @sections) ){
  $tmp_file    = "data_desc_$section.html";
  $tmp_section = "$section\_tmp.html";
  $file_name   = "classification.html";
  $subdir      = $subdirs{$file_name};

  print localtime() . "\t# Start variant classes ...\n";
  `cp $input_dir/$subdir/$file_name $tmp_file`;
  $content_before = get_content($section,'start');
  $content_after  = get_content($section,'end');
  `perl $dirname/generate_classes_table.pl -v $version -o $tmp_section -host $host -port $port -ohost $ohost -species $species`;
  $new_content = `cat $tmp_section`;
  `rm -f $tmp_section`;
  print_into_tmp_file($tmp_file,$content_before,$new_content,$content_after);

  $copy2subdir = ($no_subdir) ? '' : $subdirs{$file_name};
  copy_updated_file($copy2subdir,$file_name,$tmp_file);

  print localtime() . "\t\t> Variant classes - finished\n";
}

#### Generates the "Populations" table documentation
$section = 'populations';

if ( $sections eq "all" || grep(/^${section}$/, @sections) ){  
  $tmp_file    = "data_desc_$section.html";
  $tmp_section = "$section\_tmp.html";
  $file_name   = "populations.html";
  $subdir      = $subdirs{$file_name};

  print localtime() . "\t# Start populations ...\n";
  `cp $input_dir/$subdir/$file_name $tmp_file`;
  $content_before = get_content($section,'start');
  $content_after  = get_content($section,'end');
  `perl $dirname/generate_population_table.pl -v $version -o $tmp_section -hlist $hlist -user $user`;
  $new_content = `cat $tmp_section`;
  `rm -f $tmp_section`;
  print_into_tmp_file($tmp_file,$content_before,$new_content,$content_after);

  $copy2subdir = ($no_subdir) ? '' : $subdirs{$file_name};
  copy_updated_file($copy2subdir,$file_name,$tmp_file);

  print localtime() . "\t\t> Populations - finished\n";
}

#### Generates the "Variation sets" table documentation
$section = 'sets';

if ( $sections eq "all" || grep(/^${section}$/, @sections) ){
  $tmp_file    = "data_desc_$section.html";
  $tmp_section = "$section\_tmp.html";
  $file_name   = "sets.html";
  $subdir      = $subdirs{$file_name};

  print localtime() . "\t# Start variant sets ...\n";
  `cp $input_dir/$subdir/$file_name $tmp_file`;
  $content_before = get_content($section,'start');
  $content_after  = get_content($section,'end');
  `perl $dirname/generate_variation_set_table.pl -v $version -o $tmp_section -hlist $hlist -user $user`;
  $new_content = `cat $tmp_section`;
  `rm -f $tmp_section`;
  print_into_tmp_file($tmp_file,$content_before,$new_content,$content_after);

  $copy2subdir = ($no_subdir) ? '' : $subdirs{$file_name};
  copy_updated_file($copy2subdir,$file_name,$tmp_file);

  print localtime() . "\t\t> Variant sets - finished\n";
}

#### Generates the "Clinical significance" tables documentation
$section = 'clin_significance';

if ( $sections eq "all" || grep(/^${section}$/, @sections) ){
  $tmp_file    = "data_desc_$section.html";
  $tmp_section = "$section\_tmp.html";
  $file_name   = "phenotype_annotation.html";
  $subdir      = $subdirs{$file_name};

  print localtime() . "\t# Start clinical significance ...\n";
  `cp $input_dir/$subdir/$file_name $tmp_file`;
  $content_before = get_content($section,'start');
  $content_after  = get_content($section,'end');
  `perl $dirname/generate_clin_significance_tables.pl -v $version -o $tmp_section -host $host -port $port -species $species`;
  $new_content = `cat $tmp_section`;
  `rm -f $tmp_section`;
  print_into_tmp_file($tmp_file,$content_before,$new_content,$content_after);

  $copy2subdir = ($no_subdir) ? '' : $subdirs{$file_name};
  copy_updated_file($copy2subdir,$file_name,$tmp_file);

  print localtime() . "\t\t> Clinical significance - finished\n";
}

#### Generates the "Phenotype class" table documentation
$section = 'pheno_class';

if ( $sections eq "all" || grep(/^${section}$/, @sections) ){
  $tmp_file    = "data_desc_$section.html";
  $tmp_section = "$section\_tmp.html";
  $file_name   = "phenotype_annotation.html";
  $subdir      = $subdirs{$file_name};

  print localtime() . "\t# Start phenotype class ...\n";
  # Use output file from the clinical significance update
  if (-d "$output_dir/$subdir") {
    `cp $output_dir/$subdir/$file_name $tmp_file`;
  }
  elsif (-f "$output_dir/$file_name") {
    `cp $output_dir/$file_name $tmp_file`;
  }
  else {
    `cp $input_dir/$subdir/$file_name $tmp_file`;
  }
  $content_before = get_content($section,'start');
  $content_after  = get_content($section,'end');
  `perl $dirname/generate_pheno_class_table.pl -v $version -o $tmp_section -host $host -port $port -species $species`;
  $new_content = `cat $tmp_section`;
  `rm -f $tmp_section`;
  print_into_tmp_file($tmp_file,$content_before,$new_content,$content_after);

  $copy2subdir = ($no_subdir) ? '' : $subdirs{$file_name};
  copy_updated_file($copy2subdir,$file_name,$tmp_file);

  print localtime() . "\t\t> Phenotype class - finished\n";
}

#### Update the ontology versions
$section = "ontology";

if ( $sections eq "all" || grep(/^${section}$/, @sections) ){
  $tmp_file  = "data_desc_ontology.html";
  $file_name = "phenotype_annotation.html";
  $subdir    = $subdirs{$file_name};

  print localtime() . "\t# Start phenotype ontology ...\n";
  # Use output file from the clinical significance update
  if (-d "$output_dir/$subdir") {
    `cp $output_dir/$subdir/$file_name $tmp_file`;
  }
  elsif (-f "$output_dir/$file_name") {
    `cp $output_dir/$file_name $tmp_file`;
  }
  else {
    `cp $input_dir/$subdir/$file_name $tmp_file`;
  }
  my $sql_onto = qq{SELECT data_version FROM ontology WHERE name=? LIMIT 1};
  my $tmp_file_content = `cat $tmp_file`;
  foreach my $onto (@ontologies) {
    my $sth = get_connection_and_query("ensembl_ontology_$version", $ohost, $sql_onto, [$onto]);
    my $o_version = ($sth->fetchrow_array)[0];
    $o_version =~ s/releases\///;
    if ($o_version =~ /^(\d+):(\d+):(\d+)\s/) {
      if ($1 =~ /^\d{4}$/) {
        $o_version = "$3/$2/$1";
      }
      elsif ($3 =~ /^\d{4}$/) {
         $o_version = "$1/$2/$3";
      }
    }
    
    $tmp_file_content =~ s/<span id="$onto\_version">(\d+(-|\.)?)+<\/span>/<span id="$onto\_version">$o_version<\/span>/i;  
  }
  print_into_tmp_file($tmp_file,$tmp_file_content,'','');

  $copy2subdir = ($no_subdir) ? '' : $subdirs{$file_name};
  copy_updated_file($copy2subdir,$file_name,$tmp_file);

  print localtime() . "\t\t> Phenotype ontology - finished\n";
}

#### Create new Sources list documentation
$section = "source";

if ( $sections eq "all" || grep(/^${section}$/, @sections) ){
  $file_name = "sources_documentation.html";
  $tmp_file  = $file_name;

  print localtime() . "\t# Start sources list ...\n";
  `perl $dirname/sources2html.pl -v $version -o $tmp_file -hlist $hlist -phost $phost -config $vcf_config_file -d_dir $data_dir`;

  $copy2subdir = ($no_subdir) ? '' : $subdirs{$file_name};
  copy_updated_file($copy2subdir,$file_name,$tmp_file);

  print localtime() . "\t\t> Sources list - finished\n";
}

#### Create new Phenotype sources list documentation
$section = "pheno_source";

if ( $sections eq "all" || grep(/^${section}$/, @sections) ){  
  $file_name = "sources_phenotype_documentation.html";
  $tmp_file  = $file_name;

  print "# Phenotype sources list ...\n";
  `perl $dirname/phensources2html.pl -v $version -o $tmp_file -hlist $hlist -phost $phost`;

  $copy2subdir = ($no_subdir) ? '' : $subdirs{$file_name};
  copy_updated_file($copy2subdir,$file_name,$tmp_file);

  print localtime() . "\t\t> Phenotype sources list - finished\n";
}

#### Create the detailed species data count documentation (including SIFT and PolyPhen-2)
$section = "count";

if ( $sections eq "all" || grep(/^${section}$/, @sections) ){
  $file_name = "species_detailed_counts.html";
  $tmp_file  = $file_name;
  print localtime() . "\t# Detailed species data count ...\n";
  `perl $dirname/species_list_detailed_counts.pl -v $version -o $tmp_file -hlist $hlist -phost $phost --user ensro -p_data $p_data -config $vcf_config_file -d_dir $data_dir`;

  $copy2subdir = ($no_subdir) ? '' : $subdirs{$file_name};
  copy_updated_file($copy2subdir,$file_name,$tmp_file);

  print localtime() . "\t\t> Detailed species data count - finished\n";
}

print "\n" . localtime() . "\t>>> End of script\n";




#---------#
# Methods #
#---------#

sub get_content {
  my $section = shift;
  my $type    = shift;
  
  my $anchor = "<!-- Data $section - $type -->";
  
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

sub print_into_tmp_file {
  my $tmp    = shift;
  my $before = shift;
  my $new    = shift;
  my $after  = shift;

  open  TMP, "\t> $tmp" or die $!;
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
  my $dsn = "DBI:mysql:$dbname:$host";
  my $dbh = DBI->connect($dsn, $user, '') or die "Connection failed";

  my $sth = $dbh->prepare($sql);
  if ($params) {
    $sth->execute(join(',',@$params));
  }
  else {
    $sth->execute;
  }
  return $sth;
}

# Copy the updated temporary file to output location
sub copy_updated_file {
  my $sdir    = shift;
  my $fname   = shift;
  my $tmpfile = shift;

  my $file_path = (-d "$output_dir/$sdir") ? "$output_dir/$sdir/$fname" : "$output_dir/$fname";
  `cp $tmp_file $file_path`;
}


sub usage {
  my $msg = shift;
  print qq{
  $msg
  Usage: perl update_data_description.pl [OPTION]
  
  Update the Variation web documentation pages (under public-plugins/ensembl/htdocs/info/genome/variation/).
  
  Options:

    -help           Print this message
      
    -v              Ensembl version, e.g. 65 (Required)
    -i              Path to the input directory (Required)
    -o              Path to the output directory (Required)
    -host           Host of the human database (Required)
    -port           MySQL port of the human database (Required)
    -species        Species name. 'Homo_sapiens' by default (optional)
    -hlist          The list of host names where the new databases are stored, separated by a coma,
                    e.g. ensembldb.ensembl.org1:1234, ensembldb.ensembl.org2:1234 (Required)
    -phost          Host name where the previous databases are stored, e.g. ensembldb.ensembl.org  (Required)
    -ohost          Host name where the ontology database is stored, with the port, e.g. ensembldb.ensembl.org:1234 (Required)
    -user           MySQL user name (Required)
    -no_subdir      Doesn't copy the files onto their subdirectories    
    -dump           An output data dump file name for vcf-only species (json format) (Required)  
    -config         The location of the vcf_config.json file. By default it uses the existing one in the same repository.
    -d_dir          The directory location of where the local vcf files are stored. By default it looks in - 
                    /nfs/production/flicek/ensembl/production/ensemblftp/data_files/vertebrates 
    -p_data         Location of the json file that contain result of the previous release for comparison.
    -sections       Comma-separated list of sections of this script to run. If this options is not provided, all sections are run.
  } . "\n";
  exit(0);
}
