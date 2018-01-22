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

use strict;


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk.org>.

=cut

use warnings;

use Net::FTP;
use Getopt::Long;

my $wiggle  = 0.1; # this parameter defines the percent we allow a new file to be smaller than an old file before we complain
my $server  = 'ftp.ensembl.org';
my $release;
my $path;
my $output_dir;
my $help;

GetOptions(
    "output_dir|o=s"    => \$output_dir,
    "server=s"          => \$server,
    "release|r=s"       => \$release,       
    "path|p=s"          => \$path,
    "wiggle=s"          => \$wiggle,
    "help|h"            => \$help,
);

unless ($output_dir) {  
    warn "Must supply --output_dir argument\n" unless $help;
    $help = 1;
}

unless (defined $release) {  
    warn "Must supply --release argument\n" unless $help;
    $help = 1;
}

if ($help) {
    print "Usage: $0 --output_dir <name> --server <name> --release <num> --path <path> --wiggle <num> --help\n";
    exit(0);
}

$release--;

$path = "/pub/release-${release}/variation/gvf/" unless $path;

print "Looking for GVFs from previous release in: ${server}${path}\n";

my $ftp = Net::FTP->new($server) or die "Failed to connect to '$server:' $@";

$ftp->login;

$ftp->cwd($path) or die "Failed to cwd to '$path', is the path correct?";

my $ok_count = 0;
my $file_count = 0;

for my $species ($ftp->ls) {
    
    $ftp->cwd("$path/$species") or die "Failed to cwd to $species directory?";
    
    for my $gvf ($ftp->ls) {
    
        next unless $gvf =~ /\.gvf\.gz$/;

        $file_count++;

        my $curr_file = "$output_dir/$species/$gvf";

        if (-e $curr_file) {

            my $ftp_size = $ftp->size("$path/$species/$gvf");

            my $curr_size = -s $curr_file;

            if ($curr_size < ($ftp_size * (1 - $wiggle)) ) {
                print "$curr_file is smaller than ftp version: $curr_size vs. $ftp_size\n";
            }
            else {
                $ok_count++;
            }
        }
        else {
            print "$curr_file not found ?\n";
        }
    }
}

$ftp->quit;

print "$ok_count/$file_count GVF files look OK\n";

