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
use Test::Exception;
use Test::Warnings qw(warning :no_end_test);
use JSON;

use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Variation::DBSQL::AnnotationFileAdaptor;
our $verbose = 0;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');

my $vdb = $multi->get_DBAdaptor('variation');

# set the VCFCollection config
my $dir = $multi->curr_dir();

ok($vdb->annotation_config_file($dir.'/annotation_config.json') eq $dir.'/annotation_config.json', "DBAdaptor annotation_config_file");
my $annotation_file_adaptor = $vdb->get_AnnotationFileAdaptor;

ok($annotation_file_adaptor && $annotation_file_adaptor->isa('Bio::EnsEMBL::Variation::DBSQL::AnnotationFileAdaptor'), "isa AnnotationFileAdaptor");

my $cadd_file;
warning { $cadd_file = $annotation_file_adaptor->fetch_by_annotation_type('cadd'); };

# now we need to set the filename_template
my $temp = $cadd_file->filename_template();
$temp =~ s/###t\-root###/$dir/;
$cadd_file->filename_template($temp);
ok($cadd_file->filename_template =~ /^$dir/, "update filename_template");

my $va = $vdb->get_VariationAdaptor();
my $vfa = $vdb->get_VariationFeatureAdaptor();
#2 45417423 45417423 rs76641827 T/G
my $variation = $va->fetch_by_name('rs76641827');
my $vf = $variation->get_all_VariationFeatures()->[0];
my $cadd_scores = $cadd_file->get_scores_by_VariationFeature($vf);
ok($cadd_scores->{A} == 24.36, "CADD score for rs76641827 variant allele A");
ok($cadd_scores->{C} ==  4.36, "CADD score for rs76641827 variant allele C");
ok($cadd_scores->{G} == 14.36, "CADD score for rs76641827 variant allele G");

my $vf_cadd_scores;
warning { $vf_cadd_scores = $vf->get_cadd_scores($cadd_file->filename_template); };

ok($vf_cadd_scores->{A} == 24.36, "CADD score for VF rs76641827 variant allele A");
ok($vf_cadd_scores->{C} ==  4.36, "CADD score for VF rs76641827 variant allele C");
ok($vf_cadd_scores->{G} == 14.36, "CADD score for VF rs76641827 variant allele G");


my $gerp_file;
warning { $gerp_file = $annotation_file_adaptor->fetch_by_annotation_type('gerp'); };

# now we need to set the filename_template
$temp = $gerp_file->filename_template();
$temp =~ s/###t\-root###/$dir/;
$gerp_file->filename_template($temp);
ok($gerp_file->filename_template =~ /^$dir/, "update filename_template");

my $gerp_score = $gerp_file->get_score_by_VariationFeature($vf);
ok($gerp_score == -1, "GERP score for rs76641827");

my $vf_gerp_score;
warning { $vf_gerp_score = $vf->get_gerp_score($gerp_file->filename_template); };
ok($vf_gerp_score == -1, "GERP score for VF rs76641827");



done_testing();
