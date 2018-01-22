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

