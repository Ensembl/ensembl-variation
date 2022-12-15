#!/usr/bin/env perl
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



=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut


# Script to update the documentation page "Data description".

use strict;
use warnings;
use Getopt::Long;
use JSON;

# Print the usage instructions if run without parameters
usage() unless (scalar(@ARGV));

my ($version,$git_dir,$input_file,$output_file,$help);

GetOptions(
  'version=i' => \$version,
  'git_dir=s' => \$git_dir,
  'i=s'       => \$input_file,
  'o=s'       => \$output_file,
  'help!'     => \$help
);


usage("Ensembl release version must be specified (-version)") unless ($version);
usage("Path to the VEP_plugins directory must be specified (-git_dir)") unless ($git_dir);
usage("VEP plugin documentation input file must be specified (-i)") unless ($input_file);
usage("Output file must be specified (-o)") unless ($output_file);
usage() if ($help);

my ($content_before, $new_content, $content_after);
my $tmp_file    = 'vep_plugins_tmp.html';
my $vep_plugin_url = "https://github.com/Ensembl/VEP_plugins/";
my $vep_plugin_url_version = "https://github.com/Ensembl/VEP_plugins/blob/release/[[SPECIESDEFS::ENSEMBL_VERSION]]";
my $plugin_conf_file = "$git_dir/plugin_config.txt";
my $cpanm_url = 'https://metacpan.org/pod/';

my @files;
my %data = ();
my %data_section = ();

my %plugins_to_skip = (
  'CCDSFilter.pm' => 1,
  'CSN.pm' => 1,
  'DAS.pm' => 1,
  'GXA.pm' => 1,
  'HGVSReferenceBase.pm' => 1,
  'miRNA.pm' => 1,
  'NonSynonymousFilter.pm' => 1,
  'PolyPhen_SIFT.pm' => 1,
  'RankFilter.pm' => 1,
  'RefSeqHGVS.pm' => 1,
  'ExAC.pm' => 1,
);


my %class_colour = (
  'Conservation'              => '#02599C',
  'External ID'               => '#333333',
  'Frequency data'            => '#FF7F50',
  'HGVS'                      => '#9400D3',
  'Look up'                   => '#006400',
  'Motif'                     => '#DAA520', # goldenrod
  'Nearby features'           => '#E75480', # Dark pink
  'ND'                        => '#A9A9A9', # darkgray
  'Pathogenicity predictions' => '#1E90FF',
  'Phenotype data'            => '#22949B',
  'Publication'               => '#6A5ACD', # slateblue
  'Sequence'                  => '#5F81A9',
  'Splicing predictions'      => '#FF0000', # red
  'Structural variant data'   => '#601212',
  'Variant data'              => '#B22222',
  'Visualisation'             => '#008000' # green
);
my %class_colour_hexa = map { $_ => 1 } values %class_colour;

my @hexa_range = ('0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F');


my %plugin_extended_names = (
  'CADD'    => 'Combined Annotation Dependent Depletion',
  'CSN'     => 'Clinical Sequencing Nomenclature',
  'DAS'     => 'Distributed Annotation System',
  'G2P'     => 'gene2phenotype',
  'GO'      => 'Gene Ontology',
  'LD'      => 'Linkage Disequilibrium',
  'LoFtool' => 'Loss-of-function',
  'LOVD'    => 'Leiden Open Variation Database',
  'MPC'     => 'missense deleteriousness metric',
  'MTR'     => 'Missense Tolerance Ratio',
);

# Get information from some of the plugins
open P, "< $plugin_conf_file" or die $!;
my $plugin_name = '';
while(<P>) {
  chomp $_;
  if ($_ =~ /"key"\s=>\s"(.+)"/) {
    $plugin_name = $1;
    $data_section{$plugin_name} = ();
  }
  if ($_ =~ / "section" => "(.+)",/) {
    $data_section{$plugin_name}{'section'} = $1;
  }
}
close(P);

# Parse each plugin to extract some information
my $dh;
opendir($dh,$git_dir) or die $!;
while (my $file = readdir($dh)) {
  next if ($file !~ /\.pm$/ || $plugins_to_skip{$file});
  push (@files, $file);
  read_plugin_file($file);
}


# Print output HTML
open OUT, "> $output_file" or die $!;

my $table_content .= qq{   
  <table class="ss">
    <thead>
      <tr>
        <th>Plugin</th>
        <th style="max-width:50%">Description</th>
        <th>Category</th>
        <th>External libraries</th>
        <th>Developer</th>
      </tr>
    </thead>
    <tbody>
};
my @sorted_files = sort { lc $a cmp lc $b } @files;
my $tr_class = 'bg1';
my %plugin_class_list;

# 1 Plugin file <=> 1 row in the output table
foreach my $file (@sorted_files) {
  print "Plugin: $file\n";
  my $plugin_name  = $data{$file}{'name'};
  my $plugin_class = ($data_section{$plugin_name} && $data_section{$plugin_name}{'section'}) ? $data_section{$plugin_name}{'section'} : 'ND'; 
  my $plugin_class_colour = $class_colour{$plugin_class} ? $class_colour{$plugin_class} : get_random_colour($plugin_class);
  
  my $plugin_id = lc($plugin_class);
     $plugin_id =~ s/ /_/g;
  
  $plugin_class_list{$plugin_class} = $plugin_id;

  $table_content .= sprintf(
    '<tr class="%s plugin_row" data-category="%s">'.
    '<td><div style="font-weight:bold"><a rel="external" href="%s">%s</a></div>%s</td>'.
    '<td>%s</td>'.
    '<td><div class="vdoc_dtype_count" style="float:left;padding:2px 6px;cursor:default;background-color:%s">%s</div></td>'.
    '<td>%s</td>'.
    '<td>%s</td>'.
    '</tr>',
    $tr_class,
    $plugin_id,
    "$vep_plugin_url_version/$file",
    $plugin_name,
    ($plugin_extended_names{$plugin_name}) ? '<div style="margin-top:6px"><small>'.$plugin_extended_names{$plugin_name}.'</small></div>' : '',
    $data{$file}{'desc'},
    $plugin_class_colour,
    $plugin_class,
    ($data{$file}{'libs'} && scalar(keys(%{$data{$file}{'libs'}})) != 0) ? ((scalar(keys(%{$data{$file}{'libs'}})) > 1) ? '<ul style="padding-left:1em"><li>'.join('</li><li>',values(%{$data{$file}{'libs'}})).'</li></ul>' : (values(%{$data{$file}{'libs'}}))[0]) : '-',
    scalar(@{$data{$file}{'developer'}}) > 1 ? '<ul style="padding-left:1em"><li>'.join("</li><li>",@{$data{$file}{'developer'}})."</li></ul>" : $data{$file}{'developer'}->[0]
  );
  $tr_class = ($tr_class eq 'bg1') ? 'bg2' : 'bg1';
  $table_content .= "\n";
}
$table_content .= qq{</tbody></table>};
close(OUT);

# Add a dropdown to select the plugins to display by category
$new_content .= qq{
  <div style="margin:5px 0px 10px">
    <span style="padding-right:8px">Select categories:</span>
    <select id="select_category" onchange="javascritpt:show_hide_category()">
       <option value="all">All categories</option>};
foreach my $p_class (sort(keys(%plugin_class_list))) {
  $new_content .= sprintf(
    '<option value="%s">%s</option>',
    $plugin_class_list{$p_class},
    $p_class
  );
}
$new_content .= qq{\n    </select>\n  </div>};

$new_content .= $table_content;

# Replace the old HTML table by the new one, within the existing vep_plugins.html file
my $section = 'Plugins list';
`cp $input_file $tmp_file`;
$content_before = get_content($section,'start');
$content_after  = get_content($section,'end');
print_into_tmp_file($tmp_file,$content_before,$new_content,$content_after);

`cp $tmp_file $output_file`;
`rm -f $tmp_file`;


# Read and parse some information from the given plugin module
# Mainly parse the perldoc content
sub read_plugin_file {
  my $file = shift;
  
  my @developer = ();
  my $name = '';
  my $desc = '';
  my %libs;
  
  open F, "< $git_dir/$file" or die $!;
  while(<F>) {
    my $line = $_;
    chomp($line);

    # Get the developer information (i.e. if this has been developped by Ensembl or not)
    if ($line =~ /=head1 LICENSE/) {  
      my $license_flag = 1;
      while ($license_flag != 0) {
        $line = <F>;
        chomp($line);
        if ($line =~ /^\s*=head1/ || $line =~ /^\s*=cut/) {
          $license_flag = 0; 
        }
        else {
          @developer = ("Ensembl") if ($line =~ /Copyright\s+.+\s+EMBL-European\s+Bioinformatics\s+Institute/);
        }
      }
    }
    
    # Contact/developper information
    if ($line =~ /=head1 CONTACT/ && scalar(@developer) == 0) {
      my $contact_flag = 1;
      while ($contact_flag != 0) {
        $line = <F>;
        chomp($line);
        if ($line =~ /^\s*=head1/ || $line =~ /^\s*=cut/) {
          $contact_flag = 0; 
        }
        else {
          $line =~ s/^\s+//;
          push(@developer, $1) if ($line =~ /(.+)\s+</);
        }
      }
    }
    
    # Plugin full name (might differ from the plugin file name)
    if ($line =~ /=head1 NAME/) {
      my $name_flag = 1;
      while ($name_flag != 0) {
        $line = <F>;
        chomp($line);
        if ($line =~ /^\s*=head1/ || $line =~ /^\s*=cut/) {
          $name_flag = 0; 
        }
        else {
          $line =~ s/^\s+//;
          $name .= $1 if ($line =~ /(\S+)\s*/);
        }
      }
    }
    
    # Get the plugin description
    if ($line =~ /=head1 DESCRIPTION/) {
      my $desc_flag = 1;
      while ($desc_flag != 0) {
        $line = <F>;
        
        if ($line =~ /^\s*=head1/ || $line =~ /^\s*=cut/) {
          $desc_flag = 0;
        }
        else {
          if ($desc ne '' || $line !~ /^\s+$/) {
            if ($line =~ /^\s*\.\/(filter|vep)/ || $line =~ /^\s*--plugin/) {
              $line = '<pre class="code sh_sh">'.$line.'</pre>';
            }
            else {
              $line =~ s/^\s+// if $line =~ /\S+/;
            }
            $desc .= $line;
          }
        }
        chomp($line);
      }
    }
    
    # Get the non Ensembl Perl module dependencies
    if ($line =~ /^use\s+(.+);/) {
      my $lib = $1;
      if ($lib !~ /^(strict|warnings)/ && $lib !~ /Bio\:\:EnsEMBL\:\:/) {
        $lib =~ /^(\S+)(.*)$/;
        
        $libs{$lib} = qq{<a href="$cpanm_url$1" rel="external">$1</a>$2};
      }
    }  
  }
  close(F);
  
  # Postprocess the description content (reformatting)
  $desc = "<p>$desc</p>";
  $desc =~ s/\n\s*\n/<\/p><p>/g;
  $desc =~ s/<\/p><p><\/p>/<\/p>/g;
  $desc =~ s/\n/<br \/>/g;
  $desc =~ s/<br \/><\/p>/<\/p>/g;
  if ($desc =~ /<\/p><p>/) {
    my $lc_name = lc($name);
       $lc_name =~ s/ /_/g;
    my $desc_link = qq{ ... <a class="button" href="#$lc_name" style="padding:3px 8px 0px 8px !important;height:18px" onclick="show_hide('$lc_name');" id="a_$lc_name">more</a>};
    $desc =~ s/<\/p><p>/$desc_link<\/p><div id="div_$lc_name" style="display:none;"><p>/;
    $desc .= '</div>';
  }
  
  $data{$file} = {'name' => $name, 'desc' => $desc, 'developer' => \@developer, 'libs' => \%libs};
}


# Get the content of the existing vep_plugins.html file excluding the HTML table listing the plugins
sub get_content {
  my $section = shift;
  my $type    = shift;
  
  my $anchor = "<!-- $section - $type -->";
  
  my $line = `grep -m1 -n '$anchor' $tmp_file`;
  die "Can't find the anchor '$anchor' in the file" if (!$line || $line eq '');
  $line =~ /^(\d+):/;
  my $line_number = $1;
  my $content;
  # Beginning to start anchor
  if ($type eq 'start') {
    $content = `head -n$line_number $tmp_file`;
    if ($content !~ /$anchor(\n?)$/) {
      $content = (split("$anchor", $content))[0].$anchor;
    }
  }
  # End anchor to end
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

# Recompose a new vep_plugins.html file with the new list of plugins
sub print_into_tmp_file {
  my $tmp    = shift;
  my $before = shift;
  my $new    = shift;
  my $after  = shift;

  open  TMP, "> $tmp" or die $!;
  print TMP  $before;
  print TMP  $new;
  print TMP  $after;
}

# Generate colours for categories missing an entry in %class_colour
sub get_random_colour {
  my $plugin_class = shift;

  my $hexa = '';
  # Trying to generate darker colours (with at least 2 components starting by a number)
  while ($hexa !~ /^#\d\w\d/ && $hexa !~ /^#\w{2}\d\w\d/ && $hexa !~ /^#\d\w{3}\d/ && !$class_colour_hexa{$hexa}) {
    my @hexa_parts=('','','');
    foreach my $c (@hexa_parts){
      my $r=int(rand(16));
      $c=$hexa_range[$r].$hexa_range[$r];
    }
    $hexa = '#'.join("",@hexa_parts);
  }
  $class_colour{$plugin_class} = $hexa;
  $class_colour_hexa{$hexa} = 1;

  return $hexa;
}

sub usage {
  my $msg = shift;
  print qq{
  $msg
  Usage: perl update_web_vep_plugins_documentation.pl [OPTION]
  
  Update the page "vep_plugins.html" (under public-plugins/docs/htdocs/info/docs/tools/vep/script/).
  
  Options:

    -help           Print this message
      
    -v              Ensembl version, e.g. 65 (Required)
    -i              Input file, e.g. vep_plugins.html (Required)
    -o              Output file (Required)
    -git_dir        Path to the VEP_plugins repository - used to fetch information about the VEP plugins (Required)
  } . "\n";
  exit(0);
}
