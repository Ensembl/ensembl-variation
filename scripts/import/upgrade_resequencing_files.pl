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


use strict;
use warnings;

use Tie::File;
use Getopt::Long;

my $new_version;
my $path;

#
# This is how I run it for version 44
# bsub -q normal -W2:00 perl upgrade_resequencing_files.pl -path [path_dir] -version 44
#

GetOptions('path=s'     => \$path,
	   'version=i'  => \$new_version
	   );

usage('-version must contain a valid number') if (!$new_version);
usage('-path must contain path to directory with resequencing files files') if (!$path);


my @files = glob("$path*");
my @line;
my $call;
my $new_file;
my $previous_version = $new_version - 1;

foreach my $file (@files){
    next if($file =~ /README/);
    print "Processing.... ",$file,"\n";
    #first, unzip the file
    $call = 'gunzip '. $file;
    system($call);
    $file =~ s/\.gz//;
    #then, replace the line
    tie @line,"Tie::File",$file || die "Could not tie the $file with File::Tie: $!\n";
    $line[2] = $line[2] . " $new_version";
    #change the file name
    $new_file = $file;
    $new_file =~ s/$previous_version/$new_version/;    
    $call = 'mv ' . $file . " $new_file";
    system($call);
    #finally, compress the file again
    $call = 'gzip ' . $new_file;
    system($call);
}

sub usage {
  my $msg = shift;

  print STDERR <<EOF;

usage: perl upgrade_resequencing_files <options>

options:
    -path    <pathname>       path to the directory containing the resequencing files
    -version <new_version>    version to upgrade the files to
EOF

  die("\n$msg\n\n");
}
