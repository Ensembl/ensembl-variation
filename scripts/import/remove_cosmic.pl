
# Deletes all COSMIC data from a variation database

use strict;
use warnings;

use Getopt::Long;

use Bio::EnsEMBL::Registry;

my $registry_file;
my $verbose;
my $help;

GetOptions(
    "registry|r=s"  => \$registry_file,
    "verbose|v"     => \$verbose,
    "help|h"        => \$help,
);

unless ($registry_file) {
    print "Must supply a registry file...\n" unless $help;
    $help = 1;
}

if ($help) {
    print "Usage: $0 --registry <reg_file> --verbose --help\n";
    print "\nDeletes all COSMIC data from an Ensembl variation database\n";
    exit(0);
}

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_all($registry_file);

# COSMIC is only for human, so just get the human variation 
# database handle from the registry

my $dbh = $registry->get_adaptor(
    'human', 'variation', 'variation'
)->dbc->db_handle;

# first get our COSMIC source_id

my $get_source_sth = $dbh->prepare(qq{
    SELECT  source_id
    FROM    source
    WHERE   name = "COSMIC"
});

$get_source_sth->execute;

my ($source_id) = $get_source_sth->fetchrow_array;

die "Didn't find COSMIC source_id, does this database have COSMIC data?"
    unless defined $source_id;

print "Found COSMIC source ID: $source_id\n" if $verbose;

# now delete stuff...

# first stuff that doesn't depend on variation

# population - sample_id

$dbh->do(qq{
    DELETE  p
    FROM    population p, sample s
    WHERE   p.sample_id = s.sample_id
    AND     s.name LIKE "COSMIC:gene%"
});

# sample - name LIKE "COSMIC:gene:%" 

$dbh->do(qq{
    DELETE FROM sample WHERE name LIKE "COSMIC:gene:%"
});

print "Deleted all samples, populations and studies\n" if $verbose;

# now all the stuff that is joined to variation

# failed_variation - variation_id

$dbh->do(qq{
    DELETE  fv
    FROM    failed_variation fv, variation v
    WHERE   fv.variation_id = v.variation_id
    AND     v.source_id = $source_id
});

# flanking_sequence - variation_id

$dbh->do(qq{
    DELETE  fs
    FROM    flanking_sequence fs, variation v
    WHERE   fs.variation_id = v.variation_id
    AND     v.source_id = $source_id
});

# variation_feature - variation_id

$dbh->do(qq{
    DELETE  vf
    FROM    variation_feature vf, variation v
    WHERE   vf.variation_id = v.variation_id
    AND     v.source_id = $source_id
});

# allele - variation_id

$dbh->do(qq{
    DELETE  a
    FROM    allele a, variation v
    WHERE   a.variation_id = v.variation_id
    AND     v.source_id = $source_id
});

# variation_annotation - variation_id

$dbh->do(qq{
    DELETE  va
    FROM    variation_annotation va, variation v
    WHERE   va.variation_id = v.variation_id
    AND     v.source_id = $source_id
});

# variation_set_variation - variation_id

$dbh->do(qq{
    DELETE  vsv
    FROM    variation_set_variation vsv, variation v
    WHERE   vsv.variation_id = v.variation_id
    AND     v.source_id = $source_id
});

# variation - source_id

$dbh->do(qq{
    DELETE FROM variation WHERE source_id = $source_id
});

print "Deleted all variation associated data\n" if $verbose;

