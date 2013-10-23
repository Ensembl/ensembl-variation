# Copyright 2013 Ensembl
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

my $rca = $vdba->get_ReadCoverageAdaptor();
my $sa  = $cdba->get_SliceAdaptor();
my $pa  = $vdba->get_PopulationAdaptor();
my $ia = $vdba->get_IndividualAdaptor();

ok($rca && $rca->isa('Bio::EnsEMBL::Variation::DBSQL::ReadCoverageAdaptor'));

# test get read coverage in a region for a certain population in in a certain level
my $slice = $sa->fetch_by_region('chromosome', '1', 1, 200_000);

my $individuals = ['VENTER'];

foreach my $name (@{$individuals}){
    my $individual = shift @{$ia->fetch_all_by_name($name)};
    is($individual->name, 'VENTER', 'Test name');
    my $coverage = $rca->fetch_all_by_Slice_Individual_depth($slice, $individual, 1);
    is(scalar @$coverage, 28, 'Test coverage'); 
    my $hash;
    foreach (@$coverage) {
        my $start = $_->start;
        my $end = $_->end;
        my $slice = $_->slice->name;
        my $level = $_->level;
        #print 'sample ', $_->sample, "\n";
        my $individual = $_->individual->name;
        $hash->{$start . '-' . $end}->{'slice'} = $slice;
        $hash->{$start . '-' . $end}->{'individual'} = $individual;
        $hash->{$start . '-' . $end}->{'level'} = $level;

    }
    is($hash->{'90292-91093'}->{'slice'}, 'chromosome:GRCh37:1:1:200000:1', 'Test slice');
    is($hash->{'166096-176305'}->{'individual'}, 'VENTER', 'Test individual');
}
my @regions = @{$rca->fetch_all_regions_covered($slice, $individuals)};
is(scalar @regions, 28, 'Number of regions');

#foreach my $rc (@{$rca->fetch_all_regions_covered($slice, $individuals)}) {
#    print "range is ", $rc->[0], '-', $rc->[1], "\n";
#}

done_testing();
