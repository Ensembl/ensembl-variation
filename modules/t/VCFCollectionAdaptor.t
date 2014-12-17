# Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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
use Data::Dumper;


use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Test::MultiTestDB;


our $verbose = 0;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');

my $vdb = $multi->get_DBAdaptor('variation');
my $cdb = $multi->get_DBAdaptor('core');

my $sa = $cdb->get_SliceAdaptor();

# set the VCFCollection config
my $dir = $multi->curr_dir();
no warnings 'once';
$Bio::EnsEMBL::Variation::DBSQL::VCFCollectionAdaptor::CONFIG_FILE = $dir.'/vcf_config.json';
my $vca = $vdb->get_VCFCollectionAdaptor();

ok($vca && $vca->isa('Bio::EnsEMBL::Variation::DBSQL::VCFCollectionAdaptor'), "isa VCFCollectionAdaptor");

# fetch all
my $collections = $vca->fetch_all();
ok($collections && scalar @$collections == 1, "fetch_all count");

# fetch by ID
my $coll = $vca->fetch_by_id('1000genomes_phase1');
ok($coll && $coll->isa('Bio::EnsEMBL::Variation::VCFCollection'), "fetch_by_id isa VCFCollection");

# now we need to set the filename_template
my $temp = $coll->filename_template();
$temp =~ s/###t\-root###/$dir/;
$coll->filename_template($temp);
ok($coll->filename_template =~ /^$dir/, "update filename_template");

# get individuals
my $inds = $coll->get_all_Individuals();
ok($inds && scalar @$inds == 1092, "get_all_Individuals count 1092");
ok($inds->[0]->name eq '1000GENOMES:phase_1:HG00096', "get_all_Individuals first name is 1000GENOMES:phase_1:HG00096");

# get populations
my $pops = $coll->get_all_Populations();
ok($pops && scalar @$pops == 19, "get_all_Populations count 19");
ok($coll->has_Population('1000GENOMES:phase_1_CEU'), "has_Population 1000GENOMES:phase_1_CEU");

# fetch genotypes by slice
# my $slice = $sa->fetch_by_region('chromosome', 2, 45401130, 45421130);

# my $gts = $coll->get_all_IndividualGenotypeFeatures_by_Slice($slice);
# ok($gts && scalar @$gts == 1092, "get_all_IndividualGenotypeFeatures_by_Slice count 1092");
# ok($gts->[0]->genotype_string eq 'T|T', "get_all_IndividualGenotypeFeatures_by_Slice first genotype T|T");

# fetch genotypes by VF
my $va = $vdb->get_VariationAdaptor();
my $v  = $va->fetch_by_name('rs7569578');
my $vf = $v->get_all_VariationFeatures->[0];

$gts = $coll->get_all_IndividualGenotypeFeatures_by_VariationFeature($vf);
ok($gts && scalar @$gts == 1092, "get_all_IndividualGenotypeFeatures_by_VariationFeature count 1092");
ok($gts->[0]->genotype_string eq 'T|T', "get_all_IndividualGenotypeFeatures_by_VariationFeature first genotype T|T");

# fetch LD genotypes by slice
# my $ld_gts = $coll->_get_all_LD_genotypes_by_Slice($slice);
# ok($ld_gts && ref($ld_gts) eq 'HASH', "_get_all_LD_genotypes_by_Slice is hash");
# ok(scalar keys %$ld_gts == 374, "_get_all_LD_genotypes_by_Slice has 374 position keys");
# ok($ld_gts->{45421006} && scalar keys %{$ld_gts->{45421006}}, "_get_all_LD_genotypes_by_Slice pos 45421006 has 1092 genotypes");
# ok($ld_gts->{45411130}->{NA20811} eq 'T|A', "_get_all_LD_genotypes_by_Slice pos 45411130 ind NA20811 has genotype T|A");

done_testing();
