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

# This script looks for any minor allele frequencies of 0.5
# and then ensures that the minor allele stored in the variation
# table is not the reference allele

use strict;
use warnings;

use Getopt::Long;

use Bio::EnsEMBL::Registry;

my $registry_file;
my $help;

GetOptions(
    "registry|r=s"  => \$registry_file,
    "help|h"        => \$help,
);

unless ($registry_file) {
    print "Must supply a registry file...\n" unless $help;
    $help = 1;
}

if ($help) {
    print "Usage: $0 --registry <reg_file>\n";
    exit(0);
}

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_all($registry_file);

my $dbh = $registry->get_adaptor(
    'human', 'variation', 'variation'
)->dbc->db_handle;

my $get_vars_sth = $dbh->prepare(qq{
    SELECT  v.variation_id, v.name, v.minor_allele
    FROM    variation v, variation_feature vf
    WHERE   v.minor_allele_freq = 0.5
    AND     v.variation_id = vf.variation_id
    AND     v.minor_allele = SUBSTR(vf.allele_string,1,1)
});

my $get_syns_sth = $dbh->prepare(qq{
    SELECT  name
    FROM    variation_synonym
    WHERE   variation_id = ?
});

my $get_maf_sth = $dbh->prepare(qq{
    SELECT  allele
    FROM    maf
    WHERE   snp_id = ?
});

my $fix_maf_sth = $dbh->prepare(qq{
    UPDATE  variation
    SET     minor_allele = ?
    WHERE   variation_id = ?
});

$get_vars_sth->execute;

my $count = 0;

sub get_new_allele {
    my ($snp_id, $old_allele) = @_;

    $get_maf_sth->execute($snp_id);
    
    my $new_allele;

    while (my ($allele) = $get_maf_sth->fetchrow_array) {
        if ($allele ne $old_allele) {
            $new_allele = $allele;
            last;
        }
    }

    return $new_allele;
}

while (my ($v_id, $name, $old_allele) = $get_vars_sth->fetchrow_array) {
    
    $name =~ s/^rs//;

    my $new_allele = get_new_allele($name, $old_allele);
    
    unless ($new_allele) {
        $get_syns_sth->execute($v_id);

        while (my ($name) = $get_syns_sth->fetchrow_array) {
            next unless $name =~ /^rs/;
            $name =~ s/^rs//;
            $new_allele = get_new_allele($name, $old_allele);
            last if $new_allele;
        }
    }

    die "Didn't find alternative minor allele for variation $v_id?" unless $new_allele;

    $count++;

    $fix_maf_sth->execute($new_allele, $v_id);
}

print "Corrected $count alleles\n";

