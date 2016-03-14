# Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Test::MultiTestDB;

our $verbose = 0;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');

my $vdba = $multi->get_DBAdaptor('variation');

my $pa = $vdba->get_PopulationAdaptor;
my $sa = $vdba->get_SampleAdaptor;
my $ia = $vdba->get_IndividualAdaptor;
my ($tests, $samples, $sample, $all);

# Test SampleAdaptor
# fetch_all_by_Population 
$tests = [
    { population => '1000GENOMES:phase_1_AFR', size => 246,},
    { population => '1000GENOMES:phase_1_TSI', size => 98,},
];

foreach my $test (@$tests) {
    my $population = $pa->fetch_by_name($test->{population});
    my $samples = $sa->fetch_all_by_Population($population);
    my $size = scalar @$samples;
    is($size, $test->{size}, "Sample count for $test->{population}");
}

# fetch_all_by_name
$tests = [
    { name => 'J. CRAIG VENTER', count => 1,},
    { name => '1000GENOMES', count => 0,},
];

foreach my $test (@$tests) {
    $samples = $sa->fetch_all_by_name($test->{name});
    my $count = scalar @$samples;
    is($count, $test->{count}, "Number of returned samples for $test->{name}");
}

# fetch_by_dbID
$sample = $sa->fetch_by_dbID(8675);
is($sample->name, 'NA19122', 'Fetch by dbID 8675');

# fetch_all_by_dbID_list -- needed by web team..
my $list = [101106];
$samples = $sa->fetch_all_by_dbID_list($list);
foreach my $sample (@$samples) {
    is($sample->name, '1000GENOMES:phase_1:HG00120', "Fetch by dbID list [101106]");
}

# store
$sample = $sa->fetch_by_dbID(8675);
delete $sample->{$_} for qw(dbID name adaptor);
$sample->name('test');
ok($sa->store($sample), "store");

($sample) = @{$sa->fetch_all_by_name('test')};
ok($sample && $sample->name eq 'test', "fetch stored");

# strains

my $strains = $sa->fetch_all_strains();
is(scalar @$strains, 3, "Number of strains");

ok($sa->get_default_strains()->[0]  eq "NA12891", "default_strains");

ok($sa->get_reference_strain_name() eq "NA18635", "reference strain");

$strains = $sa->get_display_strains;
is(scalar @$strains, 3, "Number of display strains");


# storing new samples

my $individual = Bio::EnsEMBL::Variation::Individual->new(
  -name => 'test_individual',
  -type_individual => 'outbred', 
  -adaptor => $ia, 
);

$individual = $ia->store($individual);
my $individual_dbID = $individual->dbID; 

$sample = Bio::EnsEMBL::Variation::Sample->new(
  -name => 'test_sample',
  -display => 'UNDISPLAYABLE',
  -adaptor => $sa,
  -individual => $individual,
);

$sa->store($sample);
$sample = $sa->fetch_all_by_name('test_sample')->[0];

ok($individual_dbID == $sample->individual->dbID, "retrieve individual_id from associated sample");

$sample = Bio::EnsEMBL::Variation::Sample->new(
  -name => 'test_sample2',
  -adaptor => $sa,
  -display => 'UNDISPLAYABLE',
  -individual => Bio::EnsEMBL::Variation::Individual->new(
    -name => 'test_individual2',
    -type_individual => 'outbred',
    -adaptor         => $ia,
  ),
);

$sa->store($sample);
$sample = $sa->fetch_all_by_name('test_sample2')->[0];

ok($sample->individual->name eq 'test_individual2', "retrieve individual's name from associated sample");

$sample = Bio::EnsEMBL::Variation::Sample->new(
  -name => 'test_sample3',
  -adaptor => $sa,
  -display => 'UNDISPLAYABLE',
  -individual_id => $individual_dbID,
);

$sa->store($sample);
$sample = $sa->fetch_all_by_name('test_sample3')->[0];

ok($individual_dbID == $sample->individual->dbID, "Just provide individual_id for storing a new sample object");

$sample = Bio::EnsEMBL::Variation::Sample->new(
  -name => 'test_sample4',
  -adaptor => $sa,
  -display => 'UNDISPLAYABLE',
);

$sa->store($sample);
$sample = $sa->fetch_all_by_name('test_sample4')->[0];

ok($sample->individual->name eq $sample->name, "Store a new sample object without providing information on an individual.");

my $synonym = 'synonym_NA12891'; 
$samples = $sa->fetch_by_synonym($synonym);
ok(scalar @$samples == 1, 'fetch_by_synonym');

$samples = $sa->fetch_by_synonym($synonym, 'dbSNP');
ok(scalar @$samples == 1, 'fetch_by_synonym and source');

my $synonyms = $sa->fetch_synonyms(1592);
ok(scalar @$synonyms == 1, 'fetch_synonyms');
$synonyms = $sa->fetch_synonyms(1592, 'dbSNP');
ok(scalar @$synonyms == 1, 'fetch_synonyms with source name');

# _get_name_by_dbID
my $name = $sa->_get_name_by_dbID(1592);
ok($name eq 'NA12891', '_get_name_by_dbID');

my $reference_strain = $sa->fetch_reference_strain();
ok($reference_strain->name eq 'NA18635', 'fetch_reference_strain');

done_testing();
