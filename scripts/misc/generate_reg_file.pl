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

use DBI;
use FileHandle;
use Getopt::Long;

my $config = {};
GetOptions(
  $config,
  'release=i',
  'server_file=s',
  'servers=s',
  'help!',
) or die "Error: Failed to parse command line arguments\n";

usage() if ($config->{help});

$config->{server_file} ||= '/nfs/production/panda/ensembl/variation/server_file';

die "Release number is required (-release)" unless ($config->{release});
my $release = $config->{release};


my $aliases = {
  "Homo_sapiens"              => "human",
  "Mus_musculus"              => "mouse",
  "Danio_rerio"               => "zebrafish",
  "Felis_catus"               => "cat",
  "Bos_taurus"                => "cow",
  "Canis_familiaris"          => "dog",
  "Equus_caballus"            => "horse",
  "Taeniopygia_guttata"       => "zebrafinch",
  "Tetraodon_nigroviridis"    => "tetraodon",
  "Monodelphis_domestica"     => "opossum",
  "Sus_scrofa"                => "pig",
  "Ornithorhynchus_anatinus"  => "platypus",
  "Pan_troglodytes"           => "chimp",
  "Pongo_abelii"              => "orangutan",
  "Rattus_norvegicus"         => "rat",
  "Gallus_gallus"             => "chicken",
  "Drosophila_melanogaster"   => "fly",
  "Saccharomyces_cerevisiae"  => "yeast",
  "Macaca_mulatta"            => "macaque",
  "Nomascus_leucogenys"       => "gibbon",
  "Ovis_aries"                => "sheep",
  "Meleagris_gallopavo"       => "turkey",
};

my $adaptors = {
    core => {
        adaptor => 'Bio::EnsEMBL::DBSQL::DBAdaptor',
    },
    funcgen => {
        adaptor => 'Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor',
    },
    variation => {
        adaptor => 'Bio::EnsEMBL::Variation::DBSQL::DBAdaptor',
    },
};

my $fmt = q{
%s->new(
  '-species' => '%s',
  '-group'   => '%s',
  '-port'    => %d,
  '-host'    => '%s',
  '-user'    => '%s',
  '-pass'    => '%s',
  '-dbname'  => '%s',
);
};

my $alias_fmt = q{
Bio::EnsEMBL::Registry->add_alias('%s', '%s');
};

print <<END;
#!/usr/bin/env perl
use Bio::EnsEMBL::Registry;
END
for my $group (keys %$adaptors) {
  print "use " . $adaptors->{$group}->{adaptor} . ";\n";
}
print "\n";

my $dbs;

my $db_configs = {};
my $fh = FileHandle->new($config->{server_file}, 'r');
while (<$fh>) {
  chomp;
  my ($alias, $config_line) = split/\s/;
  my @pairs = split(':', $config_line);
  foreach my $pair (@pairs) {
    my ($key, $value) = split/=/, $pair;
    $db_configs->{$alias}->{$key} = $value;
  }
}

my @servers = keys %$db_configs; # set default to all servers from the server_file
@servers = split(',', $config->{servers}) if ($config->{servers});

my $group_abbrev_mappings = {
  'v' => 'variation',
  'c' => 'core',
  'f' => 'funcgen',
};

for my $server (@servers) {
  my $db_config = $db_configs->{$server};
  my $host = $db_config->{host};
  my $password = $db_config->{password} || '';
  my $user = $db_config->{user};
  my $port = $db_config->{port};
  my @groups = map { $group_abbrev_mappings->{$_} } split('', $db_config->{groups});

  my $dbh = DBI->connect(
    "DBI:mysql:host=$host;port=$port",
    $user,
    $password,
  );

  for my $group (@groups) {

    my $pattern = '%_'.$group.'_'.$release.'_%';

    my $sth = $dbh->prepare(qq{SHOW DATABASES like '$pattern'});

    $sth->execute;

    while (my ($db) = $sth->fetchrow_array) {
      next if $db =~ /master_schema/i;
      my ($species) = $db =~ /(.+)_$group/;
      $species =~ s/^(.)/uc($1)/e;
      $dbs->{$species}->{$group}->{db} = $db;
      $dbs->{$species}->{$group}->{server} = $host;
      $dbs->{$species}->{$group}->{user} = $user;
      $dbs->{$species}->{$group}->{pass} = $password;
      $dbs->{$species}->{$group}->{port} = $port;
    }
  }
}

for my $species (keys %$dbs) {
  # we only care about species with variation databases!
  next unless $dbs->{$species}->{variation};

  for my $group (keys %{ $dbs->{$species} }) {
    printf $fmt,
    $adaptors->{$group}->{adaptor},
    $species,
    $group,
    $dbs->{$species}->{$group}->{port},
    $dbs->{$species}->{$group}->{server},
    $dbs->{$species}->{$group}->{user},
    $dbs->{$species}->{$group}->{pass} || '',
    $dbs->{$species}->{$group}->{db};
  }

  if (my $alias = $aliases->{$species}) {
    printf $alias_fmt, $species, $alias;
  }
  else {
    warn "No alias for species '$species'";
  }
}

# registry files have to return a true value nowadays

print "\n1;\n";


sub usage {

  print qq{
  Print registry file for the specified servers.
  perl print_reg_file.pl -release 88 -serers staging,v1 > ensembl.registry
  Options:
    -help           Print this message
    -server_file    Default file is  /nfs/production/panda/ensembl/variation/server_file
    -servers        Comma-separated list of server short cuts e.g. v1,v2,staging,staging37,f1,f2. If not specified all servers from the server file are used
  } . "\n";
  exit(0);
}

