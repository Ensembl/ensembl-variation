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

use strict;
use warnings;

use Test::More;
use FindBin qw($Bin);

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Test::MultiTestDB;
use FileHandle;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');
my $vdba = $multi->get_DBAdaptor('variation');

#my $reg = 'Bio::EnsEMBL::Registry';
#$reg->load_all("$Bin/testnew.ensembl.registry");
#warn "using $Bin/testnew.ensembl.registry\n";
#my $vdba = $reg->get_DBAdaptor('human', 'variation');

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
#my $variation_name = 'rs7698608';
my $variation = $va->fetch_by_name($variation_name);

my $pgts = $pgta->fetch_all_by_Variation($variation);
is(scalar @$pgts, 22, 'number of pgts for variation');

my $hash;
foreach (@$pgts) {
    push @{$hash->{$_->population->name}}, $_->allele(1) . "|" . $_->allele(2);
}

## on the fly checks via compressed_genotype_var
is(join(',', sort @{$hash->{'1000GENOMES:phase_1_GBR'}}), 'C|C', 'pgts for GBR');
is(join(',', sort @{$hash->{'1000GENOMES:phase_1_LWK'}}), 'C|C,C|T', 'pgts for LWK');

my $population = $pa->fetch_by_name('PCJB:CUDAS');
my @pgts = sort {$a->dbID() <=> $b->dbID()} @{$pgta->fetch_all_by_Population($population)};

is(scalar @pgts, 3, 'number of pgts for population');
is($pgts[0]->variation()->name(), 'rs2255888', 'variation name');
is($pgts[0]->allele(1), 'T', 'allele1');
is($pgts[0]->allele(2), 'C', 'allele2');
is($pgts[0]->frequency(), 0.352252, 'frequency');
is($pgts[0]->population()->name(), 'PCJB:CUDAS', 'population name');

# store method
my $pg_hash = {
  genotype   => ['A','T'],
  frequency  => 0.12345,
  population => $pa->fetch_by_name('1000GENOMES:phase_1_GBR'),
  subsnp     => 12345,
  count      => 12345,
  variation  => $variation
};

my $pg1 = Bio::EnsEMBL::Variation::PopulationGenotype->new_fast($pg_hash);

ok($pgta->store($pg1), "store");

$pg_hash->{population} = $pa->fetch_by_name('1000GENOMES:phase_1_LWK');
my $pg2 = Bio::EnsEMBL::Variation::PopulationGenotype->new_fast($pg_hash);

my $hash2;
%$hash2 = %$pg_hash;
$hash2->{frequency} = 1 - $hash2->{frequency};
my $pg3 = Bio::EnsEMBL::Variation::PopulationGenotype->new_fast($hash2);

ok($pgta->store_multiple([$pg2, $pg3]), "store_multiple");

my $fetched = $variation->get_all_PopulationGenotypes();

ok((grep {$_->frequency && $_->frequency == 0.12345} @$fetched), "fetch stored pg 1");
ok((grep {$_->frequency && $_->frequency == 0.12345} grep {$_->population && $_->population->name eq '1000GENOMES:phase_1_LWK'} @$fetched), "fetch stored pg 2");
ok((grep {$_->frequency && $_->frequency == 0.87655} grep {$_->population && $_->population->name eq '1000GENOMES:phase_1_LWK'} @$fetched), "fetch stored pg 3");

my $fh = FileHandle->new;
my $tmpfile = "$Bin\/$$\_population_genotype.txt";
$fh->open(">$tmpfile") or die("ERROR: Could not write to $tmpfile\n");
ok($pgta->store_to_file_handle($pg1, $fh), "store_to_file_handle");
$fh->close();

open FH, $tmpfile;
my @lines = <FH>;
close FH;

unlink($tmpfile);

ok(scalar @lines == 1 && $lines[0] =~ /12345/, "check dumped data");

done_testing();
