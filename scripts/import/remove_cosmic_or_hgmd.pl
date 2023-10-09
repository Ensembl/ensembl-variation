# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2023] EMBL-European Bioinformatics Institute
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

# Deletes all COSMIC or HGMD-PUBLIC data from a variation database

use strict;
use warnings;

use Getopt::Long;

use Bio::EnsEMBL::Registry;

use IO::Handle;
STDOUT->autoflush(1); # flush STDOUT buffer immediately

my $registry_file;
my $source_name;
my $help;

GetOptions(
    "registry|r=s"  => \$registry_file,
    "source|s=s"    => \$source_name,
    "help|h"        => \$help,
);

unless ($registry_file) {
    print "Must supply a registry file...\n" unless $help;
    $help = 1;
}

unless ($source_name) {
    print "Must supply a source name...\n" unless $help;
    $help = 1;
}

if ($help) {
    print "Usage: $0 --registry <reg_file> --source <source_name> --help\n";
    print "\nDeletes all COSMIC or HGMD-PUBLIC data from an Ensembl variation database\n";
    exit(0);
}

$source_name = 'HGMD-PUBLIC' if ($source_name =~ /hgmd/i);
$source_name = uc($source_name);

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_all($registry_file);

# COSMIC and HGMD are only for human, so just get the human variation
# database handle from the registry

my $dbh = $registry->get_adaptor(
    'human', 'variation', 'variation'
)->dbc->db_handle;

# first get our source_id

my $get_source_sth = $dbh->prepare(qq{
    SELECT  source_id
    FROM    source
    WHERE   name = "$source_name"
});

$get_source_sth->execute;

my ($source_id) = $get_source_sth->fetchrow_array;

die "Didn't find $source_name source_id, does this database have $source_name data?"
    unless defined $source_id;

print "Found $source_name source ID: $source_id\n";

# now delete stuff...

# now all the stuff that is joined to variation

# failed_variation - variation_id
$dbh->do(qq{
    DELETE  fv
    FROM    failed_variation fv, variation v
    WHERE   fv.variation_id = v.variation_id
    AND     v.source_id = $source_id
});
print "- 'failed_variation' entries deleted\n";


# variation_feature, transcript_variation and MTMP_transcript_variation
my $tv_del_sth = $dbh->prepare(qq[ SELECT vf.variation_feature_id from variation_feature vf
                                          WHERE vf.source_id = $source_id;
                                         ]);
$tv_del_sth->execute() or die "Error selecting variation feature from $source_name\n";
my $tv_to_del = $tv_del_sth->fetchall_arrayref();

# Prepare SQL statements before for loop
my $del_vf_sth = $dbh->prepare(qq[ DELETE FROM variation_feature WHERE variation_feature_id = ? ]);
my $del_sth = $dbh->prepare(qq[ DELETE FROM transcript_variation WHERE variation_feature_id = ? ]);
my $del_mtmp_sth = $dbh->prepare(qq[ DELETE FROM MTMP_transcript_variation WHERE variation_feature_id = ? ]);

foreach my $to_del (@{$tv_to_del}){
  my $vf_id_del = $to_del->[0];
  $del_vf_sth->execute($vf_id_del) or die "Error deleting variation_feature_id $vf_id_del from variation_feature\n";
  $del_sth->execute($vf_id_del) or die "Error deleting variation_feature_id = $vf_id_del from transcript_variation\n";
  $del_mtmp_sth->execute($vf_id_del) or die "Error deleting variation_feature_id = $vf_id_del from MTMP_transcript_variation\n";
}
print "- 'variation_feature', 'transcript_variation' and 'MTMP_transcript_variation' entries deleted\n";

# phenotype_feature_attrib
$dbh->do(qq{
    DELETE  pfa
    FROM    phenotype_feature pf, phenotype_feature_attrib pfa
    WHERE   pf.phenotype_feature_id=pfa.phenotype_feature_id
    AND     pf.source_id = $source_id
});
print "- 'phenotype_feature_attrib' entries deleted\n";


# phenotype_feature
$dbh->do(qq{
    DELETE FROM phenotype_feature WHERE source_id = $source_id
});
print "- 'phenotype_feature' entries deleted\n";


# variation_set_variation - variation_id
$dbh->do(qq{
    DELETE  vsv
    FROM    variation_set_variation vsv, variation_set vs
    WHERE   vsv.variation_set_id = vs.variation_set_id
    AND     vs.name LIKE "\%$source_name\%"
});
print "- 'variation_set_variation' entries deleted\n";


# variation_synonym
$dbh->do(qq{
    DELETE FROM variation_synonym WHERE source_id = $source_id
});
print "- 'variation_synonym' entries deleted\n";

# variation
$dbh->do(qq{
    DELETE FROM variation WHERE source_id = $source_id
});
print "- 'variation' entries deleted\n";


print "Deleted all $source_name variation associated data\n";
