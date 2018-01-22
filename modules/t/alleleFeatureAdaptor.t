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
my $cdba = $multi->get_DBAdaptor('core');

my $afa = $vdba->get_AlleleFeatureAdaptor;
my $sa = $vdba->get_SampleAdaptor;
my $slice_a = $cdba->get_SliceAdaptor();

my $chr = '9';
my $slice = $slice_a->fetch_by_region('chromosome', $chr, 22124503, 22126503);

my $strain_name = '1000GENOMES:phase_1:NA06984';
my $sample = $sa->fetch_all_by_name($strain_name)->[0];

throws_ok { $afa->fetch_all_by_Slice('slice'); } qr/Slice arg expected/, 'Throw on wrong argument for fetch_all_by_Slice';
throws_ok { $afa->fetch_all_by_Slice($slice, 'sample'); } qr/Sample arg expected/, 'Throw on wrong argument for fetch_all_by_Slice';
throws_ok { $afa->fetch_all_by_Slice($slice, Bio::EnsEMBL::Variation::Sample->new()); } qr/Individual arg must have defined dbID/, 'Throw on wrong argument for fetch_all_by_Slice';

my $afs = $afa->fetch_all_by_Slice($slice, $sample);
my $af = $afs->[0];

my $sources = $afa->get_all_synonym_sources($af);
ok(scalar @$sources == 0, 'get_all_synonym_sources');
throws_ok { $afa->get_all_synonym_sources('AlleleFeature'); } qr/AlleleFeature argument expected/, 'Throw on wrong argument for get_all_synonym_sources';
warns_like {
  $afa->get_all_synonym_sources(Bio::EnsEMBL::Variation::AlleleFeature->new());
} qr/Not possible to get synonym sources for the AlleleFeature/, 'Warn on wrong argument for get_all_synonym_sources';
done_testing();
