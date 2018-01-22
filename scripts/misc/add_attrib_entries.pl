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

