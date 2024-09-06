#!/usr/bin/env perl
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



=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

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
  'core=s',
  'variation=s',
  'funcgen=s',
  'otherfeatures=s',
  'dbname_var=s',
  'species=s',
  'help!',
) or die "Error: Failed to parse command line arguments\n";

usage() if ($config->{help});

$config->{server_file} ||= '/nfs/production/flicek/ensembl/variation/servers.txt';
die "Error: server file not found\n" unless -e $config->{server_file};

die "Error: --release version is required\n" unless ($config->{release});
my $release = $config->{release};

my $aliases = {
  "Bos_taurus"                   => "cow",
  "Canis_familiaris"             => "dog",
  "Canis_familiarisboxer"        => "dog_boxer",
  "Capra_hircus"                 => "goat",
  "Coturnix_japonica"            => "japanese_quail",
  "Danio_rerio"                  => "zebrafish",
  "Drosophila_melanogaster"      => "fly",
  "Equus_caballus"               => "horse",
  "Felis_catus"                  => "cat",
  "Ficedula_albicollis"          => "collared_flycatcher",
  "Gallus_gallus_gca000002315v5" => "chicken",
  "Gallus_gallus"                => "chicken_old",
  "Homo_sapiens"                 => "human",
  "Macaca_fascicularis"          => "crab-eating_macaque",
  "Macaca_mulatta"               => "macaque",
  "Meleagris_gallopavo"          => "turkey",
  "Microtus_ochrogaster"         => "prairie_vole",
  "Monodelphis_domestica"        => "opossum",
  "Mus_musculus"                 => "mouse",
  "Nomascus_leucogenys"          => "gibbon",
  "Ornithorhynchus_anatinus"     => "platypus",
  "Ovis_aries_rambouillet"       => "sheep",
  "Ovis_aries"                   => "sheep_texel",
  "Pan_troglodytes"              => "chimp",
  "Pongo_abelii"                 => "orangutan",
  "Rattus_norvegicus"            => "rat",
  "Saccharomyces_cerevisiae"     => "yeast",
  "Sander_lucioperca"            => "pike-perch",
  "Sus_scrofa"                   => "pig",
  "Taeniopygia_guttata"          => "zebrafinch",
  "Tetraodon_nigroviridis"       => "tetraodon",
};

my $adaptors = {
    core => {
        adaptor => 'Bio::EnsEMBL::DBSQL::DBAdaptor',
    },
    funcgen => {
        adaptor => 'Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor',
    },
    otherfeatures => {
        adaptor => 'Bio::EnsEMBL::DBSQL::DBAdaptor',
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
Bio::EnsEMBL::Registry->add_alias('%s', '%s'); };

print <<END;
#!/usr/bin/env perl
use Bio::EnsEMBL::Registry;
END
print "use " . $adaptors->{$_}->{adaptor} . ";\n" for sort keys %$adaptors;
print "\n";

sub readServerConfig {
  my $server_file = shift;
  my $servers = {};
  my $fh = FileHandle->new($server_file, 'r');
  while (<$fh>) {
    chomp;
    next if ($_ =~ /^\s*$|^\s*#/); # skip comments or empty lines
    my ($alias, $host, $port, $user, $password) = split/\t/;
    $servers->{$alias}->{host} = $host;
    $servers->{$alias}->{port} = $port;
    $servers->{$alias}->{user} = $user;
    $servers->{$alias}->{password} = $password || '';
  }
  close $fh;
  return $servers;
}
my $servers = readServerConfig $config->{server_file};

sub getServers {
  my ($arg, @default) = @_;
  return $arg ? split ',', $arg : @default;
}
my @variation     = getServers $config->{variation},     ('v1-w', 'v3-w');
my @core          = getServers $config->{core},          ('v1',   'v3');
my @otherfeatures = getServers $config->{otherfeatures}, ('v1',   'v3');
my @funcgen       = getServers $config->{funcgen},       ('v1',   'v3');

my $dbs;
my $custom_db_exists = 0;
for my $server (keys $servers) {
  my $db_config = $servers->{$server};
  my $host = $db_config->{host};
  my $password = $db_config->{password} || '';
  my $user = $db_config->{user};
  my $port = $db_config->{port};

  # check which databases require this server
  my @groups;
  push @groups, 'variation'     if grep(/^$server$/, @variation);
  push @groups, 'core'          if grep(/^$server$/, @core);
  push @groups, 'otherfeatures' if grep(/^$server$/, @otherfeatures);
  push @groups, 'funcgen'       if grep(/^$server$/, @funcgen);
  next if !@groups;

  my $dbh = DBI->connect("DBI:mysql:host=$host;port=$port", $user, $password, );
  for my $group (@groups) {
    my $pattern = '%_'.$group.'_'.$release.'_%';
    my $sth = $dbh->prepare(qq{SHOW DATABASES like '$pattern'});
    $sth->execute;

    # Check if custom variation database exists
    my $dbname_var;
    if ($config->{dbname_var} and $group eq 'variation') {
      $dbname_var = $config->{dbname_var};
      my $sth2 = $dbh->prepare(qq{SHOW DATABASES like '$dbname_var'});
      $sth2->execute;
      $sth2->fetchrow_array and $custom_db_exists = 1;
    }
  
    while (my ($db) = $sth->fetchrow_array) {
      next if $db =~ /master_schema/i || $db =~ "${group}_.+_.+_.+";
      my ($species) = $db =~ /(.+)_$group/;
      $species =~ s/^(.)/uc($1)/e;
      $dbs->{$species}->{$group}->{db} = $dbname_var || $db;
      $dbs->{$species}->{$group}->{host} = $host;
      $dbs->{$species}->{$group}->{user} = $user;
      $dbs->{$species}->{$group}->{pass} = $password;
      $dbs->{$species}->{$group}->{port} = $port;
    }
  }
}

if ($config->{dbname_var} and !$custom_db_exists) {
  # Raise error if custom database was not found in any of the databases
  die "Database $config->{dbname_var} was not found in ",
      join(", ", @variation), "\n";
}

my @all_species;
if ($config->{species}) {
  @all_species = split ',', $config->{species};
} else {
  # if no species is user-defined, default to all species in aliases
  @all_species = keys %$aliases;
}

for my $species (@all_species) {
  # convert to scientific name if using common name
  for (keys $aliases) {
    $species = $_ if lc($species) eq $aliases->{$_};
  }
  $species = ucfirst lc $species;

  # we only care about species with variation databases!
  next unless $dbs->{$species}->{variation};

  for my $group (keys %{ $dbs->{$species} }) {
    printf $fmt,
           $adaptors->{$group}->{adaptor},
           $species,
           $group,
           $dbs->{$species}->{$group}->{port},
           $dbs->{$species}->{$group}->{host},
           $dbs->{$species}->{$group}->{user},
           $dbs->{$species}->{$group}->{pass} || '',
           $dbs->{$species}->{$group}->{db};
  }

  if (my $alias = $aliases->{$species}) {
    printf $alias_fmt, $species, $alias;
  } else {
    print "\n# No alias for species $species\n";
  }
}

print "\n\n1;\n"; # registry files have to return a true value

sub usage {

  print qq{
  Print Ensembl registry for the specified species and servers.

    # registry with write access to all variation databases
    perl generate_reg_file.pl --species mouse,human --release 108

    # registry with write access to custom variation database
    perl generate_reg_file.pl --species mouse --release 108 --dbname_var user_mus_musculus_variation_108_3_test

    # registry with read-only access to databases in v1 and v3
    perl generate_reg_file.pl --species mouse,human --release 108 -v v1,v3

    # registry with core database from staging 1, 2 and 3
    perl generate_reg_file.pl --species mouse,human --release 108 -c st1,st2,st3

  Options:

    --release       Release version (required)
    --species       Comma-separated list of species (supports scientific and
                    common names). Default: all species defined in script

    --variation     Servers for vatiation databases. Default: v1-w,v3-w
    --core          Servers for core databases. Default: v1,v3
    --funcgen       Servers for funcgen databases. Default: v1,v3
    --otherfeatures Servers for otherfeatures databases. Default: v1,v3

    --dbname_var    Custom database name for variation. Default:
                    \$species_variation_\$release_\$assembly

    --server_file   File with server details. Default:
                    /nfs/production/flicek/ensembl/variation/servers.txt
    --help          Print this message
  } . "\n";
  exit(0);
}

