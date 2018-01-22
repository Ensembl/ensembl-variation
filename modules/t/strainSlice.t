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

use Bio::EnsEMBL::Test::MultiTestDB;

use_ok('Bio::EnsEMBL::Variation::StrainSlice');

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');
my $cdb = $multi->get_DBAdaptor('core');
my $sa = $cdb->get_SliceAdaptor();
my $slice = $sa->fetch_by_region('chromosome','9',22124500,22126505);

my $strain_name = '1000GENOMES:phase_1:NA06984';
my $strain_name_2 = '1000GENOMES:phase_1:NA06986';

my $strain_slice = Bio::EnsEMBL::Variation::StrainSlice->new(
  -START   => $slice->{'start'},
  -END     => $slice->{'end'},
  -STRAND  => $slice->{'strand'},
  -ADAPTOR => $slice->adaptor(),
  -SEQ     => $slice->{'seq'},
  -SEQ_REGION_NAME => $slice->{'seq_region_name'},
  -SEQ_REGION_LENGTH => $slice->{'seq_region_length'},
  -COORD_SYSTEM    => $slice->{'coord_system'},
  -STRAIN_NAME     => $strain_name);

my $strain_slice_2 = Bio::EnsEMBL::Variation::StrainSlice->new(
  -START   => $slice->{'start'},
  -END     => $slice->{'end'},
  -STRAND  => $slice->{'strand'},
  -ADAPTOR => $slice->adaptor(),
  -SEQ     => $slice->{'seq'},
  -SEQ_REGION_NAME => $slice->{'seq_region_name'},
  -SEQ_REGION_LENGTH => $slice->{'seq_region_length'},
  -COORD_SYSTEM    => $slice->{'coord_system'},
  -STRAIN_NAME     => $strain_name_2);

my $afs = $strain_slice->get_all_AlleleFeatures;
ok(scalar @$afs == 27, 'Count AlleleFeatures for strain_1');
my $afs_2 = $strain_slice_2->get_all_AlleleFeatures;
ok(scalar @$afs_2 == 27, 'Count AlleleFeatures for strain_2');
my $diffs = $strain_slice->get_all_differences_StrainSlice($strain_slice_2);
ok(scalar @$diffs == 5, 'Differences between strain_1 and strain_2');

my $with_coverage = 1;
my $afs_with_coverage = $strain_slice->get_all_AlleleFeatures($with_coverage);
ok(scalar @$afs_with_coverage == 27, 'AlleleFeatures with coverage');

my $test_strain_name = $strain_slice->strain_name;
ok($test_strain_name eq '1000GENOMES:phase_1:NA06984', 'strain name');

my $seq = $strain_slice->seq;
ok(length($seq) == 2006, 'strain slice seq length');

my $ref_subseq  = $strain_slice->ref_subseq(1, 100, 1);
ok(length($ref_subseq) == 100, 'ref_subseq length');

my $sub_slice = $strain_slice->sub_Slice(1, 1000, 1);
ok($sub_slice && $sub_slice->isa('Bio::EnsEMBL::Variation::StrainSlice'), 'is StrainSlice');

my $afs_sub_slice = $sub_slice->get_all_AlleleFeatures();
ok(scalar @$afs_sub_slice == 14, 'AlleleFeatures for sub_slice');

my $expanded_length = $strain_slice->expanded_length;
ok($expanded_length == 2006, 'expanded_length');

my $vfs = $strain_slice->get_all_VariationFeatures();
ok(scalar @$vfs == 27, 'Count VariationFeatures');

done_testing();
