# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2017] EMBL-European Bioinformatics Institute
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
use Data::Dumper;

our $verbose = 0;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');

my $vdb = $multi->get_DBAdaptor('variation');

my $svpfa = $vdb->get_StructuralVariationPopulationFrequencyAdaptor();
my $sva   = $vdb->get_StructuralVariationAdaptor();
my $sna   = $vdb->get_SampleAdaptor();

ok($svpfa && $svpfa->isa('Bio::EnsEMBL::Variation::DBSQL::StructuralVariationPopulationFrequencyAdaptor'), "isa sv population frequency adaptor");

my $sv_name = 'esv3817090';

my $sv = $sva->fetch_by_name($sv_name);


# test fetch all by StructuralVariation

throws_ok { $svpfa->fetch_all_by_StructuralVariation('structural_variation'); } qr/Bio::EnsEMBL::Variation::StructuralVariation arg expected/, 'Throw on wrong argument for fetch_all_by_StructuralVariation';

throws_ok { $svpfa->fetch_all_by_StructuralVariation(Bio::EnsEMBL::Variation::StructuralVariation->new(-name => 'esv1')); } qr/StructuralVariation arg must have defined dbID/, 'Throw on wrong argument for fetch_all_by_StructuralVariation';

my $svpfs = $svpfa->fetch_all_by_StructuralVariation($sv);

my $sample_class_svpfs = {
  'copy_number_loss' => {
                         '19339' => 'heterozygous',
                         '19757' => 'heterozygous',
                         '19007' => 'heterozygous',
                         '20105' => 'heterozygous',
                         '19980' => 'heterozygous',
                         '19298' => 'heterozygous',
                         '19334' => 'heterozygous',
                         '18708' => 'heterozygous',
                         '19875' => 'heterozygous',
                         '19226' => 'heterozygous',
                         '19090' => 'heterozygous',
                         '18932' => 'heterozygous',
                         '19160' => 'heterozygous',
                         '19131' => 'heterozygous',
                         '19100' => 'heterozygous',
                         '20032' => 'heterozygous',
                         '19876' => 'heterozygous',
                         '19289' => 'heterozygous',
                         '19049' => 'heterozygous',
                         '18293' => 'heterozygous',
                         '19179' => 'heterozygous',
                         '19171' => 'heterozygous',
                         '20284' => 'heterozygous',
                         '20033' => 'heterozygous',
                         '19161' => 'heterozygous',
                         '19154' => 'heterozygous',
                         '19294' => 'heterozygous',
                         '20160' => 'heterozygous',
                         '18960' => 'heterozygous'
                        },
  'copy_number_gain' => {
                         '18085' => 'heterozygous',
                         '20129' => 'heterozygous',
                         '18050' => 'heterozygous'
                        }
};

is_deeply(
  $svpfs->[0]->{'samples_class'}, 
  $sample_class_svpfs,
  'Compare the "samples_class" of the first SVPF object returned'
);
ok(scalar(@$svpfs) == 14, 'Count populations');


done_testing();
