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

use Test::Exception;
use Test::More;
use FindBin qw($Bin);

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::Population;
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Test::MultiTestDB;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');
my $vdba = $multi->get_DBAdaptor('variation');

my $va = $vdba->get_VariationAdaptor;
my $pa = $vdba->get_PopulationAdaptor;
my $ia = $vdba->get_IndividualAdaptor;
my $sa = $vdba->get_SampleAdaptor;

my ($populations, $population, $all);

# fetch_all_1KG_Populations
$populations = $pa->fetch_all_1KG_Populations;
is(scalar @$populations, 54, "Number of 1000 genomes populations");

# fetch_all_HapMap_Populations
$populations = $pa->fetch_all_HapMap_Populations;
is(scalar @$populations, 12, "Number of HapMap populations");

# fetch_all_LD_Populations
$populations = $pa->fetch_all_LD_Populations;
is(scalar @$populations, 44, "Number of LD populations");

# fetch_default_LDPopulation
$population = $pa->fetch_default_LDPopulation;
is($population->name, 'CSHL-HAPMAP:HapMap-CEU', "Name for default LD population");

# fetch_by_dbID
$population = $pa->fetch_by_dbID(649);
is($population->name, 'PERLEGEN:AFD_EUR_PANEL', "Fetch by dbID 649");

$population = $pa->fetch_by_dbID(102186);
is($population->display_group_name, '1000 Genomes Project Phase 1', "Display group name for dbID fetch");


# fetch_all_by_dbID_list
my $list = [101082, 101083];
$populations = $pa->fetch_all_by_dbID_list($list);
$all = join(',', map{$_->name} sort {$a->name cmp $b->name} @$populations);
is($all, '1000GENOMES:phase_1_AFR,1000GENOMES:phase_1_EUR', "Fetch by list");
throws_ok { $pa->fetch_all_by_dbID_list('list'); } qr/list reference argument is required/, 'Throw on wrong argument for fetch_all_by_dbID_list';

# fetch_all_by_Individual
my $individual_id = 101495;
my $individual = $ia->fetch_by_dbID($individual_id);
$populations = $pa->fetch_all_by_Individual($individual);
is(scalar @$populations, 3, "Number of populations for individual HG01625");
$populations = $pa->fetch_all_by_Individual_list([$individual]);
is(scalar @$populations, 3, "Number of populations for individual HG01625");

throws_ok { $pa->fetch_all_by_Individual('individual'); } qr/Individual arg expected/, 'Throw on wrong argument for fetch_all_by_Individual';
warns_like {
  $pa->fetch_all_by_Individual(Bio::EnsEMBL::Variation::Individual->new(-name => 'individual_name'));
} qr/Individual does not have dbID, cannot retrieve Individuals/, 'Throw on wrong argument for fetch_all_by_Individual';

throws_ok { $pa->fetch_all_by_Individual_list('individual'); } qr/Listref of Bio::EnsEMBL::Variation::Individual arg expected/, 'Throw on wrong argument for fetch_all_by_Individual_list';
warns_like {
  $pa->fetch_all_by_Individual_list([Bio::EnsEMBL::Variation::Individual->new(-name => 'individual_name')]);
} qr/First Individual does not have dbID, cannot retrieve Populations/, 'Throw on wrong argument for fetch_all_by_Individual_list';

# fetch_all_by_Sample
# 1000GENOMES:phase_1:HG01625 101495 -> 1000GENOMES:phase_1_ALL,1000GENOMES:phase_1_EUR,1000GENOMES:phase_1_IBS

my @sample_ids = (101495, 101096);
my $sample = $sa->fetch_by_dbID($sample_ids[0]);

$populations = $pa->fetch_all_by_Sample($sample);

is(scalar @$populations, 3, "Number of populations for sample HG01625");
is($populations->[0]->display_group_name, '1000 Genomes Project Phase 1', "Display group name for sample HG01625's population");

throws_ok { $pa->fetch_all_by_Sample('sample'); } qr/Sample arg expected/, 'Throw on wrong argument for fetch_all_by_Sample';
warns_like {
  $pa->fetch_all_by_Sample(Bio::EnsEMBL::Variation::Sample->new(-name => 'sample_name'));
} qr/Sample does not have dbID, cannot retrieve Populations/, 'Throw on wrong argument for fetch_all_by_Sample';

throws_ok { $pa->fetch_all_by_Sample_list('sample'); } qr/Listref of Bio::EnsEMBL::Variation::Sample arg expected/, 'Throw on wrong argument for fetch_all_by_Sample_list';
warns_like {
  $pa->fetch_all_by_Sample_list([Bio::EnsEMBL::Variation::Sample->new(-name => 'sample_name')]);
} qr/First Sample does not have a dbID, cannot retrieve Populations/, 'Throw on wrong argument for fetch_all_by_Sample_list';


# fetch_all_by_Sample_list
my $samples = [];
foreach my $dbid (@sample_ids) {
   push @$samples, $sa->fetch_by_dbID($dbid); 
}
$populations = $pa->fetch_all_by_Sample_list($samples);
is(scalar @$populations, 6, "Number of populations for ssamples HG01625 and HG00109");

# fetch_by_name
$population = $pa->fetch_by_name('1000GENOMES:phase_1_IBS');
is($population->name, '1000GENOMES:phase_1_IBS', "Fetch by name 1000GENOMES:phase_1_IBS");
is($population->display_group_name, '1000 Genomes Project Phase 1', "Display group name for 1000GENOMES:phase_1_IBS");

# fetch_all_by_name_search
$populations = $pa->fetch_all_by_name_search('1000GENOMES');
is(scalar @$populations, 54, "Number of populations for fetch by name search = 1000GENOMES");

# fetch_all_by_sub_Population
$population = $pa->fetch_by_name('1000GENOMES:phase_1_IBS');
$populations = $pa->fetch_all_by_sub_Population($population);
$all = join(',', map {$_->name} sort {$a->name cmp $b->name} @$populations);
is($all, '1000GENOMES:phase_1_EUR', "Fetch by sub-population 1000GENOMES:phase_1_IBS");

throws_ok { $pa->fetch_all_by_sub_Population('population'); } qr/Population argument expected/, 'Throw on wrong argument for fetch_all_by_sub_Population';
warns_like {
  $pa->fetch_all_by_sub_Population(Bio::EnsEMBL::Variation::Population->new(-name => 'population_name'));
} qr/Cannot retrieve super populations for population without dbID/, 'Throw on wrong argument for fetch_all_by_sub_Population';

# fetch_all_by_super_Population
$population = $pa->fetch_by_name('1000GENOMES:phase_1_EUR');
# 1000GENOMES:phase_1_CEU,1000GENOMES:phase_1_FIN,1000GENOMES:phase_1_GBR,1000GENOMES:phase_1_IBS,1000GENOMES:phase_1_TSI
$populations = $pa->fetch_all_by_super_Population($population);
is(scalar @$populations, 5, "Fetch by super-population 1000GENOMES:phase_1_EUR");

throws_ok { $pa->fetch_all_by_super_Population('population'); } qr/Population argument expected/, 'Throw on wrong argument for fetch_all_by_super_Population';
warns_like {
  $pa->fetch_all_by_super_Population(Bio::EnsEMBL::Variation::Population->new(-name => 'population_name'));
} qr/Cannot retrieve sub populations for population without dbID/, 'Throw on wrong argument for fetch_all_by_super_Population';

# fetch_population_by_synonym
$populations = $pa->fetch_population_by_synonym(1372);
$all = join(',', map {$_->name} sort {$a->name cmp $b->name} @$populations);
is($all, 'PERLEGEN:AFD_AFR_PANEL', "Fetch by synonym 1372");

$populations = $pa->fetch_population_by_synonym(1372,'dbSNP');
$all = join(',', map {$_->name} sort {$a->name cmp $b->name} @$populations);
is($all, 'PERLEGEN:AFD_AFR_PANEL', "Fetch by synonym 1372 - with source");

# fetch_synonyms
my $synonyms = $pa->fetch_synonyms(650);
is(join(',', @$synonyms), 1372, "Fetch synonyms for dbID 650");

my $synonyms2 = $pa->fetch_synonyms(650,'dbSNP');
is(join(',', @$synonyms2), 1372, "Fetch synonyms for dbID 650 - with source");

# test get_dbIDs_for_population_names
my $pop_name = '1000GENOMES:phase_1_EUR';
my $dbIDs = $pa->get_dbIDs_for_population_names([$pop_name]);
my @pop_ids = map { $_ } grep {$dbIDs->{$_} eq $pop_name} keys(%$dbIDs);
ok($pop_ids[0] == 101082, 'get_dbIDs_for_population_names');

# store
$population = $pa->fetch_by_dbID(649);
delete $population->{$_} for qw(dbID name);
$population->name('test');

ok($pa->store($population), "store");
$population = $pa->fetch_by_name("test");
ok($population && $population->name eq 'test', "fetch stored");

done_testing();
