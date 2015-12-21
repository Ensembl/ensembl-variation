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

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::MappedSliceContainer;

use_ok('Bio::EnsEMBL::Variation::DBSQL::StrainSliceAdaptor');

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');
my $vdb = $multi->get_DBAdaptor('variation');
my $cdb = $multi->get_DBAdaptor('core');

my $sa = $cdb->get_SliceAdaptor();
my $strain_slice_adaptor = $vdb->get_StrainSliceAdaptor;
my $sample_adaptor = $vdb->get_SampleAdaptor;
my $slice = $sa->fetch_by_region('chromosome', '9', 22124500, 22126505);
my $msc = Bio::EnsEMBL::MappedSliceContainer->new(-SLICE => $slice);
$msc->set_StrainSliceAdaptor($strain_slice_adaptor);
my $samples = $sample_adaptor->fetch_all_by_name('1000GENOMES:phase_1:NA06984');
my $sample = $samples->[0];
$msc->attach_StrainSlice($sample);

my $mapped_slices = $msc->get_all_MappedSlices();
ok(scalar @$mapped_slices == 1, 'count MappedSlices');

my $mapped_slice = $mapped_slices->[0];
ok($mapped_slice && $mapped_slice->isa('Bio::EnsEMBL::MappedSlice'), 'is MappedSlice');

my $mapped_slice_name = $mapped_slice->name;
ok($mapped_slice_name eq 'chromosome:GRCh37:9:22124500:22126505:1#strain_1000GENOMES:phase_1:NA06984', 'mapped slice name');

done_testing();
