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
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Test::MultiTestDB;

our $verbose = 0;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');

my $vdba = $multi->get_DBAdaptor('variation');

my $pa = $vdba->get_PopulationAdaptor;
my $sa = $vdba->get_SampleAdaptor;

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

done_testing();
