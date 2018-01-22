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
use Bio::EnsEMBL::Test::MultiTestDB;
use Test::Exception;


use_ok('Bio::EnsEMBL::Variation::DBSQL::ReadCoverageAdaptor');
use_ok('Bio::EnsEMBL::Variation::ReadCoverage');

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');
my $vdb = $multi->get_DBAdaptor('variation');
my $cdb = $multi->get_DBAdaptor('core');

my $sa = $cdb->get_SliceAdaptor();
my $slice = $sa->fetch_by_region('chromosome','9',22124500,22126505);

my $sample_adapt = $vdb->get_SampleAdaptor();
my ($sample) = @{$sample_adapt->fetch_all_by_name('1000GENOMES:phase_1:NA06984')};

# get adaptor
my $rca = $vdb->get_ReadCoverageAdaptor();
ok($rca && $rca->isa('Bio::EnsEMBL::Variation::DBSQL::ReadCoverageAdaptor'), "isa adaptor");

# get coverage levels
my $levels = $rca->get_coverage_levels;
ok(ref($levels) eq 'ARRAY' && scalar @$levels == 2, "get_coverage_levels 1");
ok((grep {$_ eq '1'} @$levels) && (grep {$_ eq '2'} @$levels), "get_coverage_levels 2");

## fetch by slice

# no slice
dies_ok { $rca->fetch_all_by_Slice_Sample_depth() } "fetch_all_by_Slice_Sample_depth - no slice";

# wrong arg
dies_ok { $rca->fetch_all_by_Slice_Sample_depth($sample) } "fetch_all_by_Slice_Sample_depth - wrong arg 1";
dies_ok { $rca->fetch_all_by_Slice_Sample_depth($slice, $slice) } "fetch_all_by_Slice_Sample_depth - wrong arg 2";


# only slice
my $cov = $rca->fetch_all_by_Slice_Sample_depth($slice);
ok($cov && ref($cov) eq 'ARRAY' && scalar @$cov == 2, "fetch_all_by_Slice_Sample_depth - no args - count");
ok($cov->[0]->isa('Bio::EnsEMBL::Variation::ReadCoverage'), "fetch_all_by_Slice_Sample_depth - no args - isa ReadCoverage");

# slice and sample
$cov = $rca->fetch_all_by_Slice_Sample_depth($slice, $sample);
ok($cov && ref($cov) eq 'ARRAY' && scalar @$cov == 2, "fetch_all_by_Slice_Sample_depth - sample");

# slice and level
$cov = $rca->fetch_all_by_Slice_Sample_depth($slice, 1);
ok($cov && ref($cov) eq 'ARRAY' && scalar @$cov == 1, "fetch_all_by_Slice_Sample_depth - level");

# slice, sample, level
$cov = $rca->fetch_all_by_Slice_Sample_depth($slice, $sample, 1);
ok($cov && ref($cov) eq 'ARRAY' && scalar @$cov == 1 && $cov->[0]->level == 1, "fetch_all_by_Slice_Sample_depth - sample, level");
ok($cov->[0]->seq_region_start eq '22124503' && $cov->[0]->seq_region_end eq '22125503', "fetch_all_by_Slice_Sample_depth - coords");

# fetch all regions covered
my $ranges = $rca->fetch_all_regions_covered($slice, [$sample->name]);
my $exp = [
  [
    22124503,
    22126503
  ]
];
is_deeply($ranges, $exp, "fetch_all_regions_covered");


done_testing();
