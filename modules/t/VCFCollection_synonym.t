# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2023] EMBL-European Bioinformatics Institute
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
use Bio::EnsEMBL::Test::MultiTestDB;

use FindBin qw($Bin);
use Cwd;

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Variation::VCFCollection;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');

my $vdb = $multi->get_DBAdaptor('variation');
my $cdb = $multi->get_DBAdaptor('core');

my $dir = $multi->curr_dir();

my $sa = $cdb->get_SliceAdaptor();
ok($sa && $sa->isa('Bio::EnsEMBL::DBSQL::SliceAdaptor'), "get SliceAdaptor");

my $va = $vdb->get_VariationAdaptor;
ok($va && $va->isa('Bio::EnsEMBL::Variation::DBSQL::VariationAdaptor'), "get VariationAdaptor");

ok($vdb->vcf_config_file($dir.'/vcf_config.json') eq $dir.'/vcf_config.json', "DBAdaptor vcf_config_file");
my $vca = $vdb->get_VCFCollectionAdaptor();

ok($vca && $vca->isa('Bio::EnsEMBL::Variation::DBSQL::VCFCollectionAdaptor'), "isa VCFCollectionAdaptor");

# fetch all
my $collections = $vca->fetch_all();

# fetch by ID
my $coll = $vca->fetch_by_id('1000genomes_phase1');
ok($coll && $coll->isa('Bio::EnsEMBL::Variation::VCFCollection'), "fetch_by_id");

# Set up the filename_template to match location and file used
# for testing
my $temp = $coll->filename_template();
$temp =~ s/###t\-root###/$dir/;
$coll->filename_template($temp);
$temp =~ s/test.vcf.gz/synonym_test.vcf.gz/;
$coll->filename_template($temp);

ok($coll->filename_template =~ /^$dir/, "update filename_template - 1");
ok($coll->filename_template =~ /synonym_test.vcf.gz$/, "update filename_template - 2");

# get samples
my $samples = $coll->get_all_Samples();
ok($samples && scalar @$samples == 3, "get_all_Samples count 3");
ok($samples->[0]->name eq '1000GENOMES:phase_1:HG00096', "get_all_Samples first name is 1000GENOMES:phase_1:HG00096");

# get populations
my $pops = $coll->get_all_Populations();
ok($pops && scalar @$pops == 3, "get_all_Populations count 3");
ok($coll->has_Population('1000GENOMES:phase_1_GBR'), "has_Population 1000GENOMES:phase_1_GBR");

# fetch genotypes by VF
my $v  = $va->fetch_by_name('rs7569578');
$multi->hide('variation', 'variation_synonym');

# Temporarily add a synonym that is located in the VCF 1 base before
# the variant of interest rs7569578
$v->add_synonym('Archive dbSNP', 'synonym_1');
$va->store_synonyms($v);
my ($vf) = @{$v->get_all_VariationFeatures};
my $gts = $coll->get_all_SampleGenotypeFeatures_by_VariationFeature($vf);
ok($gts && scalar @$gts == 3, "get_all_SampleGenotypeFeatures_by_VariationFeature with nearby synonym - count 3");

# restore the previous table
$multi->restore('variation', 'variation_synonym');
done_testing();
