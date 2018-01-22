#!/usr/bin/env perl
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



=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk.org>.

=cut


# Script to dump out a table of variation sets that can be used in the documentation

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Bio::EnsEMBL::Variation::Utils::Config qw(@OVERLAP_CONSEQUENCES);


###############
### Options ###
###############
my ($web_colour_file, $web_mapping_colour, $output_file, $help);

usage() if (!scalar(@ARGV));
 
GetOptions(
     'colour_file|c=s'  => \$web_colour_file,
     'mapping_file|m=s' => \$web_mapping_colour,
     'output_file|o=s'  => \$output_file,
     'help!'            => \$help
);

if (!$output_file) {
  print "> Error! Please give an output file name using the option '-output_file'\n";
  usage();
}

if (!$web_colour_file || !-e $web_colour_file) {
  print "> Error! Please give a valid colour file name using the option '-colour_file'\n";
  usage();
}

usage() if ($help);

my $SO_BASE_LINK = 'http://www.sequenceontology.org/miso/current_release/term';


my %cons_rows;
my %consequences;
my %consequences_rank;
my %colour = get_colours_from_web();


for my $cons_set (@OVERLAP_CONSEQUENCES) {

    my $so_term  = $cons_set->{SO_term};
    my $so_acc   = $cons_set->{SO_accession};
    my $so_label = $cons_set->{label};
    my $so_desc  = $cons_set->{description};
    my $impact   = $cons_set->{impact};
    my $rank     = $cons_set->{rank};

    $so_acc = qq{<a rel="external" href="$SO_BASE_LINK/$so_acc">$so_acc</a>};

    my $row = "$so_term|$so_desc|$so_acc|$so_label|$impact";

    $cons_rows{$row} = $rank;  
}


my $cons_table = qq{
<table id="consequence_type_table" class="ss">
  <tr>
    <th style="width:5px;padding-left:0px;padding-right:0px;text-align:center">*</th>
    <th>SO term</th>
    <th>SO description</th>
    <th>SO accession</th>
    <th>Display term</th>
    <th><span class="_ht ht" title="Classification of the level of severity of the consequence type">IMPACT</span></th>
  </tr>\n};

my $bg = '';
my $border_top = ';border-top:1px solid #FFF';
my $not_first = 0;

for my $row (sort {$cons_rows{$a} <=> $cons_rows{$b} || $a cmp $b} keys(%cons_rows)) {
  my $SO_term = (split(/\|/, $row))[0];
  $row =~ s/\|/<\/td>\n    <td>/g;
    
  # Fetch the group colour
  $row =~ /^(\S+)</;
  my $c = ($colour{lc($1)}) ? $colour{lc($1)} : $colour{'default'};
  my $border = ($not_first == 1) ? $border_top : '';
  
  my $cons_line = ($not_first == 0) ? '' : qq{  </tr>\n};
  $cons_line .= qq{  <tr$bg id="$SO_term">\n};
  # Consequence colour
  $cons_line .= (defined($c)) ? qq{    <td style="padding:0px;margin:0px;background-color:$c$border"></td>} : qq{    <td></td>};
  # Consequence data
  $cons_line .= qq{    <td>$row</td>\n};
  $not_first = 1;
  $bg = ($bg eq '') ? qq{ class="bg2"} : '';
    
  $cons_table .= qq{$cons_line  </tr>\n};
}

$cons_table .= qq{</table>\n};
$cons_table .= qq{<p><b>*</b> Corresponding colours for the Ensembl web displays.<p>\n};

open OUT, "> $output_file" or die $!;
print OUT $cons_table;
close(OUT);


# Retrieve the variation consequence colours from the COLOUR.ini file
sub get_colours_from_web {
  my %webcolours;
  my $var_flag = 0;
  open F, "$web_colour_file" or die $!;
  while(<F>) {
    chomp $_;
    if ($_ =~ /\[variation\]/) {
      $var_flag = 1;
      next;
    }
    next if ($var_flag == 0);
    
    if ($_ =~ /^(\S+)\s*=\s*(\S+)\s*/) {
      my $cons = $1;
      my $c    = (split(';',$2))[0]; 
      
      
      # Convert the colour name into hexadecimal code
      if (-e $web_mapping_colour) {
        my $line = `grep -w "'$c'" $web_mapping_colour`;
        $c = '#'.$1 if ($line =~ /=>\s+'(\S+)'/);
      }
      
      $webcolours{$cons} = $c;
    }
    
    last if ($_ eq '');
  }
  close(F);
  return (scalar(keys(%webcolours))) ? %webcolours : %colour;
}


sub usage {
  
  print qq{
  Usage: perl generate_consequence_table.pl [OPTION]
  
  Put all variation consequence information into an HTML table.
  
  Options:

    -help           Print this message
      
    -output_file    An HTML output file name (Required)      
    -colour_file    If you want to use directly the colours from the web colours configuration file
                    instead of the almost-up-to-date-colour-hash \%colour hash. (optional)
                    Usually, you can find the colour configuration file in:
                    ensembl-webcode/conf/ini-files/COLOUR.ini 
    
    -mapping_file   Web module to map the colour names to the corresponding hexadecimal code. (optional)
                    Useful because some colour names are internal to Ensembl and won't be displayed in the 
                    documentation pages (i.e. not using the perl modules).
                    The module NamedColours.pm can be find in:
                    ensembl-webcode/modules/EnsEMBL/Draw/Utils/NamedColours.pm
  } . "\n";
  exit(0);
}

