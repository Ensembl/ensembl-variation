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

use strict;
use warnings;

use Test::More;
use FindBin qw($Bin);

use Bio::EnsEMBL::Registry;

my $reg = 'Bio::EnsEMBL::Registry';
$reg->load_all("$Bin/test.ensembl.registry.72");

my $cdba = $reg->get_DBAdaptor('human', 'core');
my $vdba = $reg->get_DBAdaptor('human', 'variation');

my $pgta = $vdba->get_PopulationGenotypeAdaptor();

my $va = $vdba->get_VariationAdaptor();
my $pa = $vdba->get_PopulationAdaptor();


my $pgt = $pgta->fetch_by_dbID(12);

is($pgt->variation->name, 'rs1162', 'variation name');
is($pgt->allele(1), 'A', 'allele 1');
is($pgt->allele(2), 'A', 'allele 2');
is($pgt->frequency(), 0.4, 'frequency');
is($pgt->population->name, 'KWOK:HydatidiformMoles', 'population name');


my $variation_name = 'rs144235347';
my $variation = $va->fetch_by_name($variation_name);

my $pgts = $pgta->fetch_all_by_Variation($variation);
is(scalar @$pgts, 22, 'number of pgts for variation');

my $hash;
foreach (@$pgts) {
    push @{$hash->{$_->population->name}}, $_->allele(1) . "|" . $_->allele(2);
}

is(join(',', sort @{$hash->{'1000GENOMES:phase_1_GBR'}}), 'C|C', 'pgts for GBR');
is(join(',', sort @{$hash->{'1000GENOMES:phase_1_LWK'}}), 'C|C,C|T', 'pgts for LWK');

my $population = $pa->fetch_by_name('PCJB:CUDAS');
my @pgts = sort {$a->dbID() <=> $b->dbID()} @{$pgta->fetch_all_by_Population($population)};

is(scalar @pgts, 9, 'number of pgts for population');
is($pgts[0]->dbID(), 8670757, 'dbID');
is($pgts[0]->variation()->name(), 'rs2255888', 'variation name');
is($pgts[0]->allele(1), 'T', 'allele1');
is($pgts[0]->allele(2), 'C', 'allele2');
is($pgts[0]->frequency(), 0.352252, 'frequency');
is($pgts[0]->population()->name(), 'PCJB:CUDAS', 'population name');

done_testing();
