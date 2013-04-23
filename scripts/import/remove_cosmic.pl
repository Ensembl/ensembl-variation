

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

# population - population_id

$dbh->do(qq{  DELETE FROM population WHERE name LIKE "COSMIC:gene%"});

print "Deleted all populations\n" if $verbose;

# now all the stuff that is joined to variation

# failed_variation - variation_id

$dbh->do(qq{
    DELETE  fv
    FROM    failed_variation fv, variation v
    WHERE   fv.variation_id = v.variation_id
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

# phenotype_feature_attrib - variation_id

$dbh->do(qq{
    DELETE  pfa
    FROM    phenotype_feature pf, phenotype_feature_attrib pfa
    WHERE   pf.phenotype_feature_id=pfa.phenotype_feature_id
    AND     pf.source_id = $source_id
});

# phenotype_feature - variation_id

$dbh->do(qq{
    DELETE  pf
    FROM    phenotype_feature pf
    WHERE   pf.source_id = $source_id
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

