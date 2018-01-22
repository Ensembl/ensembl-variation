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

use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Test::MultiTestDB;

our $verbose = 0;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');

my $vdb = $multi->get_DBAdaptor('variation');


my $pfpma = $vdb->get_ProteinFunctionPredictionMatrixAdaptor();


ok($pfpma && $pfpma->isa('Bio::EnsEMBL::Variation::DBSQL::ProteinFunctionPredictionMatrixAdaptor'), "isa");


my $pfpm = $pfpma->fetch_by_analysis_translation_md5('sift', 'baf58fc2ba1c23683f1a67902f3a6989');

my @prediction = $pfpm->prediction_from_matrix(33, 'K');

ok($pfpm->analysis() eq 'sift',                                       "analysis type");
ok($pfpm->translation_md5() eq 'baf58fc2ba1c23683f1a67902f3a6989',    "translation_md5");

ok($prediction[0] eq "tolerated",                                     "sift prediction");
ok($prediction[1] eq 0.09,                                            "sift score");

ok($pfpm->evidence_for_prediction( 30, 'conservation_score') eq 2.78, "sift conservation_score");
ok($pfpm->evidence_for_prediction( 30, 'sequence_number') eq 156,     "sift sequence_number");



my $pfpm2 = $pfpma->fetch_by_analysis_translation_md5('polyphen_humvar', 'baf58fc2ba1c23683f1a67902f3a6989');

my @prediction2 = $pfpm2->prediction_from_matrix(33, 'K');

ok($prediction2[0] eq "benign",                "polyphen prediction");
ok($prediction2[1] eq 0.027,                   "polyphen score");



my $pfpm3 = $pfpma->fetch_by_analysis_translation_md5('polyphen_humdiv', 'baf58fc2ba1c23683f1a67902f3a6989');

my @prediction3 = $pfpm2->prediction_from_matrix(33, 'K');

ok($prediction3[0] eq "benign",                "polyphen humdiv prediction");
ok($prediction3[1] eq 0.027,                   "polyphen humdiv score");



my $pfpm4 = $pfpma->fetch_sift_predictions_by_translation_md5( 'baf58fc2ba1c23683f1a67902f3a6989');

my @prediction4 = $pfpm->prediction_from_matrix(33, 'K');

ok($pfpm4->analysis() eq 'sift',               "sift analysis type");
ok($prediction4[0] eq "tolerated",             "sift prediction");
ok($prediction4[1] eq 0.09,                    "sift score");



my $pfpm5 = $pfpma->fetch_polyphen_predictions_by_translation_md5('baf58fc2ba1c23683f1a67902f3a6989');

my @prediction5 = $pfpm5->prediction_from_matrix(33, 'K');

ok($pfpm5->analysis() eq 'polyphen',           "analysis type");
ok($pfpm5->sub_analysis() eq 'humvar',         "sub analysis type");
ok($prediction5[0] eq "benign",                "polyphen humvar prediction");
ok($prediction5[1] eq 0.027,                   "polyphen humvar score");



$pfpm->peptide_length(300); 
ok($pfpm->peptide_length() == 300,             "peptide_length");



done_testing();
