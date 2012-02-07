use strict;
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

unless (defined $release && !$help) {  
    warn "Must supply --release argument\n" unless $help;
    $help = 1;
}

if ($help) {
    print "Usage: $0 --output_dir <name> --server <name> --release <num> --path <path> --wiggle <%> --help\n";
    exit(0);
}

$release--;

$path = "/pub/release-${release}/variation/gvf/" unless $path;

my $ftp = Net::FTP->new($server) or die "Failed to connect to '$server:' $@";

$ftp->login;

$ftp->cwd($path) or die "Failed to cwd to '$path', is the path correct?";

for my $species ($ftp->ls) {
    
    $ftp->cwd("$path/$species") or die "Failed to cwd to $species directory?";
    
    for my $gvf ($ftp->ls) {
    
        next unless $gvf =~ /\.gvf\.gz$/;

        my $curr_file = "$output_dir/$species/$gvf";

        if (-e $curr_file) {

            my $ftp_size = $ftp->size("$path/$species/$gvf");

            my $curr_size = -s $curr_file;

            if ($curr_size < ($ftp_size * (1 - $wiggle)) ) {
                print "$curr_file is smaller than ftp version: $curr_size vs. $ftp_size\n";
            }
        }
        else {
            print "$curr_file not found ?\n";
        }
    }
}

$ftp->quit;

