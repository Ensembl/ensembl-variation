#!/use/bin/env perl

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
    "Pongo_pygmaeus"            => "orangutan",
    "Rattus_norvegicus"         => "rat",
    "Gallus_gallus"             => "chicken",
    "Drosophila_melanogaster"   => "fly",
    "Saccharomyces_cerevisiae"  => "yeast",
};

my $groups = {
    core => {
        user    => 'ensro',
        adaptor => 'Bio::EnsEMBL::DBSQL::DBAdaptor',
    },
    funcgen => {
        user    => 'ensro',
        adaptor => 'Bio::EnsEMBL::Variation::DBSQL::DBAdaptor',
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
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Registry;
END

my @species;

for my $server (@servers) {

    my $dbh = DBI->connect(
        "DBI:mysql:host=$server;port=3306",
        'ensro',
        '',
    );
    
    my $sth = $dbh->prepare(qq{
        SHOW DATABASES like '%_variation_${release}_%'
    });#'
    
    $sth->execute;

    while (my ($db) = $sth->fetchrow_array) {

        my ($species) = $db =~ /(.+)_variation/;
        
        $species =~ s/^(.)/uc($1)/e;
        
        for my $group (qw(core variation)) {
            
            my $group_db = $db;
            
            $group_db =~ s/_variation_/_${group}_/;
            
            printf $fmt,
                $groups->{$group}->{adaptor},
                $species,
                $group,
                $default_port,
                $server,
                $groups->{$group}->{user},
                $groups->{$group}->{pass} || '',
                $group_db;
        }
        
        if (my $alias = $aliases->{$species}) {
            printf $alias_fmt, $species, $alias;
        }
        else {
            warn "No alias for species '$species'";
        }
    }
}

