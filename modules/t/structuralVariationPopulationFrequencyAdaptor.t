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
                         '19339' => 1,
                         '19757' => 1,
                         '19007' => 1,
                         '20105' => 1,
                         '19980' => 1,
                         '19298' => 1,
                         '19334' => 1,
                         '18708' => 1,
                         '19875' => 1,
                         '19226' => 1,
                         '19090' => 1,
                         '18932' => 1,
                         '19160' => 1,
                         '19131' => 1,
                         '19100' => 1,
                         '20032' => 1,
                         '19876' => 1,
                         '19289' => 1,
                         '19049' => 1,
                         '18293' => 1,
                         '19179' => 1,
                         '19171' => 1,
                         '20284' => 1,
                         '20033' => 1,
                         '19161' => 1,
                         '19154' => 1,
                         '19294' => 1,
                         '20160' => 1,
                         '18960' => 1
                        },
  'copy_number_gain' => {
                         '18085' => 1,
                         '20129' => 1,
                         '18050' => 1
                        }
};

is_deeply(
  $svpfs->[0]->{'samples_class'}, 
  $sample_class_svpfs,
  'Compare the "samples_class" of the first SVPF object returned'
);
ok(scalar(@$svpfs) == 14, 'Count populations');


done_testing();
