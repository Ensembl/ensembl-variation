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

# Script to generate GFF files containing GO terms for all species


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut


use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::VEP::Config;
use DBI;
use strict;
use Getopt::Long;
use GO;

###############
### Options ###
###############
my ($output_dir,$e_version,$user,$hname);
my ($force,$quiet, $species, $help);

## EG options
my ($site, $etype);

usage() if (!scalar(@ARGV));

GetOptions(
     'output_dir|o=s'  => \$output_dir,
     'force|f'         => \$force,
     'quiet|q'         => \$quiet,
     "species|s=s"     => \$species,
     'v=s'             => \$e_version,
     'user|u=s'        => \$user,
     'hname=s'         => \$hname,
     'help!'           => \$help,
);

if ($help) {
  usage();
} elsif (!$e_version) {
  print STDERR "> Error! Please give an Ensembl version, using the option '-v' \n";
  usage();
} elsif (!$output_dir) {
  print STDERR "> Error! Please give an output directory for the dumps, using the option '-o' \n";
  usage();
} elsif (!$user) {
  print STDERR "> Error! Please give a database user name, using the option '-user' \n";
  usage();
} elsif (!$hname) {
  print STDERR "> Error! Please give the host name where the databases are stored using the option '-hname'\n";
  usage();
}

# Settings
my $database = "";

#VEP-config
my $cfg = Bio::EnsEMBL::VEP::Config->new({testing => 'hello'});
#registry data:
my ($host, $port) = split /\:/, $hname;
my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db(
    -host       => $host,
    -port       => $port,
    -user       => $user,
    -pass       => '',
    -db_version => $e_version,
 );

my %only_species = map { $_ => 1 } split(",",$species) if defined $species;

# Get variation databases 
my $sql = qq{SHOW DATABASES LIKE '%variation_$e_version%'};
my $sth_h = get_connection_and_query($database, $hname, $sql, 1);

# Loop over variation databases
my $total = 0;
while (my ($dbname) = $sth_h->fetchrow_array) {
  next if ($dbname =~ /^master_schema/ || $dbname =~ /private/);

  $dbname =~ /^(.+)_variation_.+_(.+)/;
  my $s_name = $1;
  my $assembly = $2;
  # Add GRCh prefix as expected by web VEP
  $assembly = "GRCh" . $assembly if $s_name eq 'homo_sapiens';

  # Check core database
  my $dbcore = $dbname;
  $dbcore =~ s/variation/core/g;
  next if defined $species & !$only_species{$s_name};
  print STDERR "\n# $s_name - [ $dbcore ]\n";

  # Skip Monodelphis domestica that gives an error when creating tabix index
  if ($s_name eq "monodelphis_domestica") {
    warn "Skipping $s_name because of large genome size\n";
    next;
  }

  # Check for GO terms in core database
  my $sqlGO = qq{
    SELECT * FROM xref x
    JOIN external_db db ON x.external_db_id = db.external_db_id
    WHERE db.db_name = "GO" AND x.dbprimary_acc LIKE "GO:%" };
  my $sth_go = get_connection_and_query($dbcore, $hname, $sqlGO, 1);
  my ($count) = $sth_go->fetchrow_array();
  print "Found $count GO terms\n";
  next unless $count > 0;

  # Run GO.pm plugin to export GFF-formatted GO terms
  my %reg_config = (reg => $registry, species => $s_name, assembly => $assembly,
                    quiet => $quiet );
  my @params = ($output_dir);
  my %db_config = (config => \%reg_config, params => \@params, match => 'transcript' );
  mkdir $output_dir unless -d $output_dir;
  my $dumpFile = $output_dir . "/" . GO::_prepare_filename(\%db_config, $registry);
  
  if ( -e $dumpFile & !$force ) {
    warn "Skipping $dumpFile. Use --force to override the existing file\n";
    next;
  }
  GO::_generate_gff(\%db_config, $dumpFile);
  $total++;
}
$sth_h->finish;
print "\n\n$total GFF file(s) containing GO terms were dumped to $output_dir\n";

# Connects and execute a query
sub get_connection_and_query {
  my $dbname  = shift;
  my $hname   = shift;
  my $sql     = shift;
  my $rtn_sth = shift;

  my ($host, $port) = split /\:/, $hname;

  # DBI connection 
  my $dsn = "DBI:mysql:$dbname:$host:$port";
  my $dbh = DBI->connect($dsn, $user) or die "Connection failed";

  my $sth = $dbh->prepare($sql);
  $sth->execute;

  if ($rtn_sth) {
    return $sth;
  }
  else {
    $sth->finish;
  }
}


sub usage {

  print qq{
  Usage: perl dump_GO_terms.pl [OPTION]

  Dumping GO terms from Ensembl database to GFF files.

  Options:

    -help             Print this message

    -output_dir|o     Output directory for the GFF dumps (Required)
    -force|f          Force to override any pre-existing dumps
    -quiet|q          Suppress warning messages. Not used by default

    -species|s        Only export the species in the list,
                        e.g. homo_sapiens,sus_scrofa

    -user|u           Database login user name (Required)
    -v                Ensembl version, e.g. 96 (Required)
    -hname            The host name (with port) where the databases are stored,
                        e.g. ensembldb.ensembl.org1:3334 (Required)
  } . "\n";
  exit(0);
}
