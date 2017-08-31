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

our $verbose = 0;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');

my $vdb = $multi->get_DBAdaptor('variation');

my $svpfa = $vdb->get_StructuralVariationPopulationFrequencyAdaptor();
my $sva   = $vdb->get_StructuralVariationAdaptor();
my $sna   = $vdb->get_SampleAdaptor();

ok($svpfa && $svpfa->isa('Bio::EnsEMBL::Variation::DBSQL::StructuralVariationPopulationFrequencyAdaptor'), "isa sv population frequency adaptor");

my $sample_name = 'NA18635'; 
my $ind_name    = 'NA18635'; 
my $ind_gender  = 'Male';
my $sv_name     = 'esv3817090';
my $study_name  = 'estd214';



my $sv = $sva->fetch_by_name($sv_name);


# test fetch all by StructuralVariation

throws_ok { $svpfa->fetch_all_by_StructuralVariation('structural_variation'); } qr/Bio::EnsEMBL::Variation::StructuralVariation arg expected/, 'Throw on wrong argument for fetch_all_by_StructuralVariation';
throws_ok { $svpfa->fetch_all_by_StructuralVariation(Bio::EnsEMBL::Variation::StructuralVariation->new(-name => 'esv1')); } qr/StructuralVariation arg must have defined dbID/, 'Throw on wrong argument for fetch_all_by_StructuralVariation';

my $svpfs = $svpfa->fetch_all_by_StructuralVariation($sv);

ok(scalar(@$svpfs) == 14, 'Count populations');


done_testing();
