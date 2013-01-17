use strict;

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

