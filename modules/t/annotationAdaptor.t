# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2024] EMBL-European Bioinformatics Institute
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

use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Test::MultiTestDB;
our $verbose = 0;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');

my $vdb = $multi->get_DBAdaptor('variation');

my $dir = $multi->curr_dir();

ok($vdb->vcf_config_file($dir.'/vcf_config.json') eq $dir.'/vcf_config.json', "DBAdaptor vcf_config_file");
my $gerp_annotation_adaptor = $vdb->get_GERPAnnotationAdaptor;
my $cadd_annotation_adaptor = $vdb->get_CADDAnnotationAdaptor;

ok($gerp_annotation_adaptor && $gerp_annotation_adaptor->isa('Bio::EnsEMBL::Variation::DBSQL::GERPAnnotationAdaptor'), "isa GERPAnnotationAdaptor");

my $gerp_annotation = $gerp_annotation_adaptor->fetch_by_id('70_mammals.gerp_conservation_score');
my $cadd_annotation = $cadd_annotation_adaptor->fetch_by_id('CADD');

# now we need to set the filename_template
$vdb->gerp_root_dir($dir);
ok($vdb->gerp_root_dir eq $dir, "set and get gerp_root_dir");
my $temp = $gerp_annotation->filename_template();
$temp =~ s/###t\-root###/$dir/;

$gerp_annotation->filename_template($temp);
ok($gerp_annotation->filename_template =~ /^$dir/, "update GERPAnnotation filename_template");

$temp = $cadd_annotation->filename_template();
$temp =~ s/###t\-root###/$dir/;
$cadd_annotation->filename_template($temp);
ok($cadd_annotation->filename_template =~ /^$dir/, "update CADDAnnotation filename_template");

my $va = $vdb->get_VariationAdaptor();
my $vfa = $vdb->get_VariationFeatureAdaptor();
#2 45417423 45417423 rs76641827 T/G
my $variation = $va->fetch_by_name('rs76641827');
my $vf = $variation->get_all_VariationFeatures()->[0];

my $gerp_score = $gerp_annotation->get_score_by_VariationFeature($vf);
ok($gerp_score == -1.12, "GERP score for rs76641827");

my $vf_gerp_score = $vf->get_gerp_score($gerp_annotation->filename_template);
my ($id, $score) = %{$vf_gerp_score};
ok($score == -1.12, "GERP score for VF rs76641827");
ok($id eq '70_mammals.gerp_conservation_score', "GERP annotation file name");

my $cadd_scores = $cadd_annotation->get_all_scores_by_VariationFeature($vf);
ok($cadd_scores->{A} == 24.36, "CADD score for rs76641827 variant allele A");
ok($cadd_scores->{C} ==  4.36, "CADD score for rs76641827 variant allele C");
ok($cadd_scores->{G} == 14.36, "CADD score for rs76641827 variant allele G");

my $vf_cadd_scores = $vf->get_all_cadd_scores($cadd_annotation->filename_template);
my ($cadd_annotation_id, $scores) = %{$vf_cadd_scores};
ok($cadd_annotation_id eq 'CADD', "CADD annotation file name");
ok($scores->{A} == 24.36, "CADD score for VF rs76641827 variant allele A");
ok($scores->{C} ==  4.36, "CADD score for VF rs76641827 variant allele C");
ok($scores->{G} == 14.36, "CADD score for VF rs76641827 variant allele G");

# test some corner cases
$variation = $va->fetch_by_name('rs2299222');
$vf = $variation->get_all_VariationFeatures()->[0];
$vf_gerp_score = $vf->get_gerp_score($gerp_annotation->filename_template);
($id, $score) = %{$vf_gerp_score};
ok(! defined $score, "Returns undef if GERP score not available in annotation file");

# CADD for insertions rs35107173
$variation = $va->fetch_by_name('rs35107173');
$vf = $variation->get_all_VariationFeatures()->[0];

warns_like {
  $vf_cadd_scores = $vf->get_all_cadd_scores($cadd_annotation->filename_template);
} qr/Can only calculate CADD scores for variants of length 1/, 'Warn if input variant is an insertion';

# location not in CADD file:
$variation = $va->fetch_by_name('rs199476127');
$vf = $variation->get_all_VariationFeatures()->[0];
$vf_cadd_scores = $vf->get_all_cadd_scores($cadd_annotation->filename_template);
($cadd_annotation_id, $scores) = %{$vf_cadd_scores};
ok($cadd_annotation_id eq 'CADD', "CADD annotation file name");
ok(scalar keys %$scores == 0, "Return empty if location is not in CADD file");

# GERP for insertions rs70937952
$variation = $va->fetch_by_name('rs70937952');
$vf = $variation->get_all_VariationFeatures()->[0];
$vf_gerp_score = $vf->get_gerp_score($gerp_annotation->filename_template);
($id, $score) = %{$vf_gerp_score};
ok($score == -8.37, "GERP score for VF rs70937952 insertion");

done_testing();
