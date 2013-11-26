#!/usr/bin/env perl
# Copyright [1999-2013] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk.org>.

=cut


use strict;
use warnings;

use DBI;

my $release = $ARGV[0] || die "Usage: $0 <release_num>\n";

my @servers = qw(ens-staging ens-staging2);

my $default_port = 3306;

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
};

my $groups = {
    core => {
        user    => 'ensro',
        adaptor => 'Bio::EnsEMBL::DBSQL::DBAdaptor',
    },
    funcgen => {
        user    => 'ensro',
        adaptor => 'Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor',
    },
    variation => {
        user    => 'ensadmin',
        pass    => 'ensembl',
        adaptor => 'Bio::EnsEMBL::Variation::DBSQL::DBAdaptor',
    },
};

my $fmt = q{
%s->new(
    '-species'  => '%s',
    '-group'    => '%s',
    '-port'     => %d,
    '-host'     => '%s',
    '-user'     => '%s',
    '-pass'     => '%s',
    '-dbname'   => '%s',
);
};

my $alias_fmt = q{
Bio::EnsEMBL::Registry->add_alias('%s', '%s');
};

print <<END;
#!/usr/bin/env perl
use Bio::EnsEMBL::Registry;
END

for my $group (keys %$groups) {
    print "use ".$groups->{$group}->{adaptor}.";\n";
}

print "\n";

my $dbs;

for my $server (@servers) {

    my $dbh = DBI->connect(
        "DBI:mysql:host=$server;port=3306",
        'ensro',
        '',
    );
   
    for my $group (keys %$groups) {

        my $pattern = '%_'.$group.'_'.$release.'_%';

        my $sth = $dbh->prepare(qq{SHOW DATABASES like '$pattern'});
        
        $sth->execute;

        while (my ($db) = $sth->fetchrow_array) {
            next if $db =~ /master_schema/i;
            my ($species) = $db =~ /(.+)_$group/;
            $species =~ s/^(.)/uc($1)/e;
            $dbs->{$species}->{$group}->{db} = $db;
            $dbs->{$species}->{$group}->{server} = $server;
        }
    }
}

for my $species (keys %$dbs) {

    # we only care about species with variation databases!
    next unless $dbs->{$species}->{variation};

    for my $group (keys %{ $dbs->{$species} }) {
        printf $fmt,
                $groups->{$group}->{adaptor},
                $species,
                $group,
                $default_port,
                $dbs->{$species}->{$group}->{server},
                $groups->{$group}->{user},
                $groups->{$group}->{pass} || '',
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

