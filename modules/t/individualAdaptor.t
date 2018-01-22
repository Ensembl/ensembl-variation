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
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Test::MultiTestDB;

our $verbose = 0;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');

my $vdba = $multi->get_DBAdaptor('variation');

my $pa = $vdba->get_PopulationAdaptor;
my $ia = $vdba->get_IndividualAdaptor;

my ($tests, $individuals, $individual, $all);

# Test IndividualAdaptor
# fetch_all_by_Population 
$tests = [
  { population => '1000GENOMES:phase_1_AFR', size => 246,},
  { population => '1000GENOMES:phase_1_TSI', size => 98,},
];

foreach my $test (@$tests) {
  my $population = $pa->fetch_by_name($test->{population});
  my $individuals = $ia->fetch_all_by_Population($population);
  my $size = scalar @$individuals;
  is($size, $test->{size}, "Individual count for $test->{population}");
}
throws_ok { $ia->fetch_all_by_Population('population'); } qr/Population arg expected/, 'Throw on wrong argument for fetch_all_by_Population';
warns_like { 
  $ia->fetch_all_by_Population(Bio::EnsEMBL::Variation::Population->new(-name => 'population_name'));
} qr/Population does not have dbID/, 'Throw on wrong argument for fetch_all_by_Population';
# fetch_all_by_name
$tests = [
  { name => 'J. CRAIG VENTER', count => 1,},
  { name => '1000GENOMES', count => 0,},
];

foreach my $test (@$tests) {
  $individuals = $ia->fetch_all_by_name($test->{name});
  my $count = scalar @$individuals;
  is($count, $test->{count}, "Number of returned individuals for $test->{name}");
}

$individuals = $ia->fetch_all_by_name_list(['1000GENOMES:phase_1:NA19660']);
ok($individuals->[0]->name eq '1000GENOMES:phase_1:NA19660', 'fetch_all_by_name_list');
throws_ok { $ia->fetch_all_by_name_list('list'); } qr/list reference argument is required/, 'Throw on wrong argument';

# fetch_by_dbID
$individual = $ia->fetch_by_dbID(8675);
is($individual->name, 'NA19122', 'Fetch by dbID 8675');

# fetch_all_by_dbID_list -- needed by web team..
my $list = [101106];
$individuals = $ia->fetch_all_by_dbID_list($list);
foreach my $individual (@$individuals) {
    is($individual->name, '1000GENOMES:phase_1:HG00120', "Fetch by dbID list [101106]");
}

# fetch_all_by_parent_Individual
my $parent = $ia->fetch_by_dbID(101961);
my $mother = $ia->fetch_by_dbID(101960);
is($parent->name, '1000GENOMES:phase_1:NA19661', "Parent name is 1000GENOMES:phase_1:NA19661");

my $populations = $parent->get_all_Populations();
$all = join(',', map {$_->name} sort {$a->name cmp $b->name} @$populations);
is($all,'1000GENOMES:phase_1_ALL,1000GENOMES:phase_1_AMR,1000GENOMES:phase_1_MXL', "Populations membership for 1000GENOMES:phase_1:NA19661");

my $children = $ia->fetch_all_by_parent_Individual($parent);
$all = join(',', map {$_->name} sort {$a->name cmp $b->name} @$children);
is($all, '1000GENOMES:phase_1:NA19685', "All children for 1000GENOMES:phase_1:NA19661");

$children = $ia->fetch_all_by_parent_Individual($mother);
$all = join(',', map {$_->name} sort {$a->name cmp $b->name} @$children);
is($all, '1000GENOMES:phase_1:NA19685', "All children for 1000GENOMES:phase_1:NA19660");

throws_ok { $ia->fetch_all_by_parent_Individual('individual'); } qr/Individual argument expected/, 'Throw on wrong argument for fetch_all_by_parent_Individual';
warns_like { 
  $ia->fetch_all_by_parent_Individual(Bio::EnsEMBL::Variation::Individual->new(-name => 'parent'));
} qr/Cannot fetch child Individuals for parent without dbID/, 'Throw on wrong argument for fetch_all_by_parent_Individual';

my $father_individual = Bio::EnsEMBL::Variation::Individual->new(-name => 'father', -description => 'father', -gender => 'Unknown');
my $mother_individual = Bio::EnsEMBL::Variation::Individual->new(-name => 'mother', -description => 'mother', -gender => 'Unknown');
my $father = $ia->store($father_individual);
my $father_individual_id = $father->dbID;
$mother = $ia->store($mother_individual);
my $mother_individual_id = $mother->dbID;
my $child = Bio::EnsEMBL::Variation::Individual->new(-name => 'child', -description => 'child', -father_individual_id => $father_individual_id, -mother_individual_id => $mother_individual_id,);
$child->father_Individual($father);
$child->mother_Individual($mother);
$child = $ia->store($child);
my $child_individual_id = $child->dbID;
$children = $ia->fetch_all_by_parent_Individual($father);
$all = join(',', map {$_->name} sort {$a->name cmp $b->name} @$children);
is($all, 'child', "fetch children by parent without specified gender");
my $dbh = $ia->dbc->db_handle;
$dbh->do(qq{DELETE FROM individual WHERE individual_id IN ($father_individual_id, $mother_individual_id, $child_individual_id);}) or die $dbh->errstr;

# synonyms 
#my $ind =  $ia->fetch_synonyms(101101);
ok( $ia->fetch_synonyms(101101)->[0] eq "fred", "fetch synonym");
ok( $ia->fetch_synonyms(101101, 'dbSNP')->[0] eq "fred", "fetch synonym by dbID and source");

$individuals = $ia->fetch_individual_by_synonym("fred");
ok($individuals->[0]->name eq "1000GENOMES:phase_1:HG00114", "fetch by synonym");

$individuals = $ia->fetch_individual_by_synonym("fred", 'dbSNP');
ok($individuals->[0]->name eq "1000GENOMES:phase_1:HG00114", "fetch by synonym and source");

# store
$individual = $ia->fetch_by_dbID(8675);
delete $individual->{$_} for qw(dbID name adaptor);
$individual->name('test');
ok($ia->store($individual), "store");

($individual) = @{$ia->fetch_all_by_name('test')};
ok($individual && $individual->name eq 'test', "fetch stored");

my $individual_name = $ia->_get_name_by_dbID(101960);
ok($individual_name eq '1000GENOMES:phase_1:NA19660', '_get_name_by_dbID');

done_testing();
