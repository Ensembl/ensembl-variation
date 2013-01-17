#!/usr/bin/env perl

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
