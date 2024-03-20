#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2024] EMBL-European Bioinformatics Institute
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


# Script to check or update the vep documentation pages.

use strict;
use warnings;
use DBI;
use Getopt::Long;
use File::Basename;
use Mojo::DOM;

use Bio::EnsEMBL::VEP::Config;

# Print the usage instructions if run without parameters
usage() unless (scalar(@ARGV));

my ($help,$input_dir);

GetOptions(
  'help!'       => \$help,
  'i=s'         => \$input_dir
);

usage("input directory must be specified") unless $input_dir;

my ($section,$tmp_file,$file_name,$subdir);
my %subdirs = ( 'vep_options.html'              => 'script',
              );

$file_name   = "vep_options.html";
$tmp_file    = "tmp_$file_name.html";
$subdir      = $subdirs{$file_name};

print localtime() . "\t# Check incompatibilities in vep options page...\n";
`cp $input_dir/$subdir/$file_name $tmp_file`;

my %INCOMPATIBLE = %Bio::EnsEMBL::VEP::Config::INCOMPATIBLE;
my $incompatible_matric = {};
foreach (keys %INCOMPATIBLE) {
  $incompatible_matric->{$_} = [] unless $incompatible_matric->{$_};
  push $incompatible_matric->{$_}, @{ $INCOMPATIBLE{$_} };
  foreach my $reverse_keys (@{ $INCOMPATIBLE{$_} }) {
    $incompatible_matric->{$reverse_keys} = [] unless $incompatible_matric->{$reverse_keys};
    push @{ $incompatible_matric->{$reverse_keys} }, $_;
  }
}

open my $IN, "<", $tmp_file or die "Could not open $tmp_file - $!";
{
  local $/;
  my $content = <$IN>;
  # print $content;
  my $dom = Mojo::DOM->new($content);
  use Data::Dumper;
  for (@{ $dom->find('table') }) {
    my $incompatible_idx;
    my $idx = 0;
    for (@{ $_->find('th') }) {
      $incompatible_idx = $idx if $_->content eq "Incompatible with";
      $idx++;
    }

    if (defined $incompatible_idx){
      for (@{ $_->find('tr') }) {
        my $option;
        if ($_->{'id'}) {
          $option = $_->{'id'};
          $option =~ s/opt_//;
        }

        if (defined $option) {
          my $incompatible_col = $_->find('td')->[$incompatible_idx];
          my $incompatibilities = [];
          for (@{ $incompatible_col->find('a') }){
            push @{ $incompatibilities }, substr($_->content, 2);
          }

          
          # missing incompatibilities
          for my $opt (@{ $incompatible_matric->{$option} }){
            print "$option missing incompatibility: $opt\n" unless grep(/^$opt$/, @{ $incompatibilities });
          }

          # redundant incompatibilities
          for my $opt (@{ $incompatibilities }){
            print "$option redundant incompatibility: $opt\n" unless grep(/^$opt$/, @{ $incompatible_matric->{$option} });
          }
        }
      }
    }
  }
}

sub usage {
  my $msg = shift;
  print qq{
  $msg
  Usage: perl update_web_vep_documentation.pl [OPTION]
  
  Check or update the VEP web documentation pages (under public-plugins/info/docs/tools/vep/).
  Currently only check for required updates in incompatible options in vep options page.
  
  Options:

    -help           Print this message
    -i              Path to the input directory (Required)
  } . "\n";
  exit(0);
}
