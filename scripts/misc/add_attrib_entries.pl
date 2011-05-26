use strict;
use warnings;

use DBI;
use Getopt::Long;

my $host;
my $port = 3306;
my $user;
my $pass;
my $pattern;
my $attrib_file;
my $dry_run;
my $help;

GetOptions(
    "host=s"        => \$host,
    "port=i"        => \$port,
    "user=s"        => \$user,
    "pass=s"        => \$pass,
    "pattern=s"     => \$pattern,
    "attrib_file=s" => \$attrib_file,
    "dry_run"       => \$dry_run,
    "help|h"        => \$help,
);

unless ($host && $user && $pattern && $attrib_file) {
    print "Missing required parameter...\n" unless $help;
    $help = 1;
}

if ($help) {
    print "Usage: $0 --host <host> --port <port> --user <user> --pass <pass> --pattern <pattern> --attrib_file <attrib_sql_file> --dry_run --help\n";
    exit(0);
}

unless (-e $attrib_file) {
    die "Can't see file: '$attrib_file'\n";
}

my $dbh = DBI->connect(
    "DBI:mysql:host=$host;port=$port",
    $user,
    $pass,
);

my $sth = $dbh->prepare("SHOW DATABASES LIKE '$pattern'");

$sth->execute;

while (my ($db) = $sth->fetchrow_array) {
    my $cmd = "mysql --host=$host --port=$port --user=$user --pass=$pass --database=$db < $attrib_file";
    print "CMD: $cmd\n";
    unless ($dry_run) {
        my $rc = system($cmd);
        if ($rc != 0) {
            die "command failed...";
        }
    }
}

