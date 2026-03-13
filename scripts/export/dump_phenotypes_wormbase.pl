# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2026] EMBL-European Bioinformatics Institute
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

# Script to generate an HTML page containing the variant sources of each species


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut


use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::VEP::Config;
use ProductionMysql;
use DBI;
use strict;
use Getopt::Long;
use Phenotypes;

###############
### Options ###
###############
my ($output_dir,$e_version,$p_version,$user,$hname);
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
     'pv=s'            => \$p_version,
     'user|u=s'        => \$user,
     'hname=s'         => \$hname,
     'help!'           => \$help,
);

if ($help) {
  usage();
} elsif (!$p_version) {
  print STDERR "> Error! Please give a WormBase ParaSite version, using the option '-pv' \n";
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

my $sql = qq{SHOW DATABASES LIKE '%variation_${p_version}_${e_version}%'};
my $sth_h = get_connection_and_query($database, $hname, $sql, 1);

my $sqlPF = qq{SELECT count(*) FROM phenotype_feature};

# Loop over databases
while (my ($dbname) = $sth_h->fetchrow_array) {

  next if ($dbname =~ /^master_schema/ || $dbname =~ /private/);

  next if ($dbname !~ /^[a-z]+_[a-z]+_[a-z0-9]+_variation_\d+_\d+_\d+$/i);

  my $cdbname = `sed 's/variation/core/g' <<< $dbname`;
  #print "Core: ".$cdbname."\n";

  my $ftpfilename = ProductionMysql::core_db_to_local_ftp_filename(ProductionMysql->staging, $cdbname);

  $dbname =~ /^(.+)_variation_.+_(.+)/;
  my $s_name = $1;
  my $assembly = $2;

  next if defined $species & !$only_species{$s_name};
  print STDERR "# $s_name - [ $dbname ]\n";

  ## count phenotype_feature(s)
  my $sth_pf = get_connection_and_query($dbname, $hname, $sqlPF, 1);
  my ($count) = $sth_pf->fetchrow_array();
  print "Found phenotype features: ", $count, "\n";
  next unless $count > 0;

  my $dumpFile = sprintf($output_dir."/".$ftpfilename."."."phenotypes.wb", $s_name, $e_version);

  die ("Phenotype file $dumpFile already exists. Specify a different output folder or use --force to override the existing one\n") if -e $dumpFile & !$force ;

  print "writing: $dumpFile \n" unless $quiet;

  my %reg_config = (reg => $registry, species => $s_name, quiet => $quiet, dbname => $dbname, cdbname => $cdbname);
  my %db_config = (config => \%reg_config);

  #Run Phenotypes.pm plugin to export GVF formated Phenotypes
  Phenotypes::generate_phenotype_gff(\%db_config, $dumpFile);
}
$sth_h->finish;



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
  Usage: perl dump_phenotypes_wormbase.pl [OPTION]

  Dumping to GVF of phenotype data for WormBase ParaSite species.

  Options:

    -help             Print this message

    -output_dir|o     Output directory for the GVF dumps (Required)
    -force|f          Force to override any pre-existing dumps
    -quiet|q          Suppress warning messages.Not used by default

    -species|s        Only export the species in the list,
                        e.g. caenorhabditis_elegans,schistosoma_mansoni

    -user|u           Database login user name (Required)
    -v                Ensembl version, e.g. 96 (Required)
    -pv               WormBase ParaSite Version, e.g. 16 (Required)       
    -hname            The host name (with port) where the databases are stored,
                        e.g. ensembldb.ensembl.org1:3334 (Required)
  } . "\n";
  exit(0);
}
