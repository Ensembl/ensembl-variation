# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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
use Bio::EnsEMBL::Variation::Population;
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Test::MultiTestDB;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');
my $vdba = $multi->get_DBAdaptor('variation');


my $pa = $vdba->get_PopulationAdaptor;
my $ia = $vdba->get_IndividualAdaptor;

my ($populations, $population, $all);

# fetch_all_1KG_Populations
$populations = $pa->fetch_all_1KG_Populations;
is(scalar @$populations, 22, "Number of 1000 genomes populations");

# fetch_all_HapMap_Populations
$populations = $pa->fetch_all_HapMap_Populations;
is(scalar @$populations, 12, "Number of HapMap populations");

# fetch_all_LD_Populations
$populations = $pa->fetch_all_LD_Populations;
is(scalar @$populations, 18, "Number of LD populations");

# fetch_default_LDPopulation
$population = $pa->fetch_default_LDPopulation;
is($population->name, 'CSHL-HAPMAP:HapMap-CEU', "Name for default LD population");

# fetch_by_dbID
$population = $pa->fetch_by_dbID(649);
is($population->name, 'PERLEGEN:AFD_EUR_PANEL', "Fetch by dbID 649");

# fetch_all_by_dbID_list
my $list = [101082, 101083];
$populations = $pa->fetch_all_by_dbID_list($list);
$all = join(',', map{$_->name} sort {$a->name cmp $b->name} @$populations);
is($all, '1000GENOMES:phase_1_AFR,1000GENOMES:phase_1_EUR', "Fetch by list");

# fetch_all_by_Individual
# 1000GENOMES:phase_1:HG01625 101495 -> 1000GENOMES:phase_1_ALL,1000GENOMES:phase_1_EUR,1000GENOMES:phase_1_IBS

my @individual_dbids = (101495, 101096);
my $individual = $ia->fetch_by_dbID($individual_dbids[0]);
$populations = $pa->fetch_all_by_Individual($individual);
is(scalar @$populations, 3, "Number of populations for individual HG01625");

# fetch_all_by_Individual_list
my $individuals = [];
foreach my $dbid (@individual_dbids) {
   push @$individuals, $ia->fetch_by_dbID($dbid); 
}
$populations = $pa->fetch_all_by_Individual_list($individuals);
is(scalar @$populations, 6, "Number of populations for individuals HG01625 and HG00109");

# fetch_by_name
$population = $pa->fetch_by_name('1000GENOMES:phase_1_IBS');
is($population->name, '1000GENOMES:phase_1_IBS', "Fetch by name 1000GENOMES:phase_1_IBS");

# fetch_all_by_name_search
$populations = $pa->fetch_all_by_name_search('1000GENOMES');
is(scalar @$populations, 22, "Number of populations for fetch by name search = 1000GENOMES");

# fetch_all_by_sub_Population
$population = $pa->fetch_by_name('1000GENOMES:phase_1_IBS');
$populations = $pa->fetch_all_by_sub_Population($population);
$all = join(',', map {$_->name} sort {$a->name cmp $b->name} @$populations);
is($all, '1000GENOMES:phase_1_EUR', "Fetch by sub-population 1000GENOMES:phase_1_IBS");

# fetch_all_by_super_Population
$population = $pa->fetch_by_name('1000GENOMES:phase_1_EUR');
# 1000GENOMES:phase_1_CEU,1000GENOMES:phase_1_FIN,1000GENOMES:phase_1_GBR,1000GENOMES:phase_1_IBS,1000GENOMES:phase_1_TSI
$populations = $pa->fetch_all_by_super_Population($population);
is(scalar @$populations, 5, "Fetch by super-population 1000GENOMES:phase_1_EUR");

# fetch_population_by_synonym
$populations = $pa->fetch_population_by_synonym(1372);
$all = join(',', map {$_->name} sort {$a->name cmp $b->name} @$populations);
is($all, 'PERLEGEN:AFD_AFR_PANEL', "Fetch by synonym 1372");

# fetch_synonyms
my $synonyms = $pa->fetch_synonyms(650);
is(join(',', @$synonyms), 1372, "Fetch synonyms for dbID 650");

# fetch_tag_Population
#my $v = $va->fetch_by_name('rs205621');
#my $vfs = $vfa->fetch_all_by_Variation($v);
#$populations = $pa->fetch_tag_Population($vfs->[0]);

# fetch_tagged_Population

# get_sample_id_for_population_names

# store
$population = $pa->fetch_by_dbID(649);
delete $population->{$_} for qw(dbID name);
$population->name('test');

ok($pa->store($population), "store");
$population = $pa->fetch_by_name("test");
ok($population && $population->name eq 'test', "fetch stored");

done_testing();
