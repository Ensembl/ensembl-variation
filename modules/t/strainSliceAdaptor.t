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

my $sample_name = '1000GENOMES:phase_1:NA06984';
my $mapped_slices = $strain_slice_adaptor->fetch_by_name($msc, $sample_name);
my $mapped_slice = $mapped_slices->[0];
my $pairs = $mapped_slice->get_all_Slice_Mapper_pairs();
ok(scalar @$pairs == 1, 'get_all_Slice_Mapper_pairs');
my $pair = $pairs->[0];
my ($slice_mapped, $mapper) = @$pair;

my $sgfa = $vdb->get_SampleGenotypeFeatureAdaptor;

my $sample = $sample_adaptor->fetch_all_by_name('1000GENOMES:phase_1:HG00096')->[0];
my $dir = $multi->curr_dir();
ok($vdb->vcf_config_file($dir.'/vcf_config.json') eq $dir.'/vcf_config.json', "DBAdaptor vcf_config_file");

my $vca = $vdb->get_VCFCollectionAdaptor();
my $coll = $vca->fetch_by_id('1000genomes_phase1');
# now we need to set the filename_template
my $temp = $coll->filename_template();
$temp =~ s/###t\-root###/$dir/;
$coll->filename_template($temp);

$sgfa->db->use_vcf(2);

$slice = $sa->fetch_by_region('chromosome', '2', 45401280, 45421006);

my $genotypes = $sgfa->fetch_all_by_Slice($slice, $sample);
$msc = Bio::EnsEMBL::MappedSliceContainer->new(-SLICE => $slice, -EXPANDED => 1);
$mapped_slices = $strain_slice_adaptor->fetch_by_name($msc, '1000GENOMES:phase_1:HG00096');
$mapped_slice = $mapped_slices->[0];
my $seq = $mapped_slice->seq(1);
ok( length $seq == 19727, 'mapped_slice length');

done_testing();
