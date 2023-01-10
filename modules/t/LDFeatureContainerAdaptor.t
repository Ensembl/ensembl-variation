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
use Test::Warnings qw(warning :no_end_test);
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Variation::VariationFeature;
our $verbose = 0;

use_ok('Bio::EnsEMBL::Variation::LDFeatureContainer');
use_ok('Bio::EnsEMBL::Variation::DBSQL::LDFeatureContainerAdaptor');

## examples to be updated & added to test-genome-DBs files

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();

my $vdb = $multi->get_DBAdaptor('variation');
my $db = $multi->get_DBAdaptor('core');

$vdb->dnadb($db);

my $ldfca = $vdb->get_LDFeatureContainerAdaptor();
my $ldContainer;

ok($ldfca && $ldfca->isa('Bio::EnsEMBL::Variation::DBSQL::LDFeatureContainerAdaptor'), "get adaptor");
$ldfca->db->use_vcf(1);

my $sa = $db->get_SliceAdaptor();
my $vfa = $vdb->get_VariationFeatureAdaptor();
my $va = $vdb->get_VariationAdaptor();
my $pa = $vdb->get_PopulationAdaptor();

my $dir = $multi->curr_dir();
ok($vdb->vcf_config_file($dir.'/ld_vcf_config.json') eq $dir.'/ld_vcf_config.json', "DBAdaptor vcf_config_file");

my $vca = $vdb->get_VCFCollectionAdaptor();
my $coll = $vca->fetch_by_id('ld');
# now we need to set the filename_template
my $temp = $coll->filename_template();
$temp =~ s/###t\-root###/$dir/;
$coll->filename_template($temp);

# fetch_by_Slice
my $slice = $sa->fetch_by_region('chromosome', '9', 22124503, 22126503);
my $population = $pa->fetch_by_name('1000GENOMES:phase_1_ASW');
my $population_no_genotypes_in_vcf = $pa->fetch_by_name('PGA-UW-FHCRC:HSP_GENO_PANEL');

my $ldfc;
warning { $ldfc = $ldfca->fetch_by_Slice($slice, $population_no_genotypes_in_vcf); };
my $ld_values = $ldfc->get_all_ld_values;
cmp_ok(scalar @$ld_values, '==', 0, "Return empty container if population is not present in any VCF file");

$ldfc = $ldfca->fetch_by_Slice($slice, $population);
$ld_values = $ldfc->get_all_ld_values;

my $results_fetch_by_Slice = [
{ population_id => 102179, vf1_name => 'rs79944118', vf2_name => 'rs1333049', r2 => 0.087731, d_prime => 0.999965, test_name => 'ld_values for rs79944118 and rs1333049' },
{ population_id => 102179, vf1_name => 'rs10757279', vf2_name => 'rs1333048', r2 => 0.630300, d_prime => 0.999998, test_name => 'ld_values for rs10757279 and rs1333048' },
{ population_id => 102179, vf1_name => 'rs10757279', vf2_name => 'rs1333050', r2 => 0.402075, d_prime => 0.662126, test_name => 'ld_values for rs10757279 and rs1333050' },
{ population_id => 102179, vf1_name => 'rs1333047', vf2_name => 'rs1333049', r2 => 0.050071, d_prime => 0.999871, test_name => 'ld_values for rs1333047 and rs1333049' },
{ population_id => 102179, vf1_name => 'rs1333048', vf2_name => 'rs1333050', r2 => 0.233557, d_prime => 0.635636, test_name => 'ld_values for rs1333048 and rs1333050' },
{ population_id => 102179, vf1_name => 'rs1333048', vf2_name => 'rs1333049', r2 => 0.684916, d_prime => 0.999999, test_name => 'ld_values for rs1333048 and rs1333049' },
{ population_id => 102179, vf1_name => 'rs10757279', vf2_name => 'rs1333049', r2 => 0.614303, d_prime => 0.817026, test_name => 'ld_values for rs10757279 and rs1333049' },
{ population_id => 102179, vf1_name => 'rs79944118', vf2_name => 'rs1333048', r2 => 0.060082, d_prime => 0.999914, test_name => 'ld_values for rs79944118 and rs1333048' },
{ population_id => 102179, vf1_name => 'rs4977575', vf2_name => 'rs1333049', r2 => 0.050071, d_prime => 0.999871, test_name => 'ld_values for rs4977575 and rs1333049' },
{ population_id => 102179, vf1_name => 'rs4977575', vf2_name => 'rs72655407', r2 => 0.063754, d_prime => 0.999996, test_name => 'ld_values for rs4977575 and rs72655407' },
{ population_id => 102179, vf1_name => 'rs79944118', vf2_name => 'rs186966498', r2 => 0.108295, d_prime => 0.46932, test_name => 'ld_values for rs79944118 and rs186966498' },
{ population_id => 102179, vf1_name => 'rs1333047', vf2_name => 'rs72655407', r2 => 0.063754, d_prime => 0.999996, test_name => 'ld_values for rs1333047 and rs72655407' },
{ population_id => 102179, vf1_name => 'rs1333049', vf2_name => 'rs1333050', r2 => 0.354847, d_prime => 0.648413, test_name => 'ld_values for rs1333049 and rs1333050' },
{ population_id => 102179, vf1_name => 'rs1333047', vf2_name => 'rs4977575', r2 => 1, d_prime => 1, test_name => 'ld_values for rs1333047 and rs4977575' },
];

foreach my $result (@$results_fetch_by_Slice) {
  my $vf1 = ($va->fetch_by_name($result->{vf1_name})->get_all_VariationFeatures)->[0];
  my $vf2 = ($va->fetch_by_name($result->{vf2_name})->get_all_VariationFeatures)->[0];
  my $population_id = $result->{population_id};
  my $r2 = $result->{r2};
  my $d_prime = $result->{d_prime};
  my $test_name = $result->{test_name};  
  cmp_ok($ldfc->get_r_square($vf1, $vf2, $population_id), '==', $r2, "$test_name r2");
  cmp_ok($ldfc->get_d_prime($vf1, $vf2, $population_id), '==', $d_prime, "$test_name d_prime");
}

#fetch_all_by_Variation
my $v = $va->fetch_by_name('rs1333049');
my $ldfcs = $ldfca->fetch_all_by_Variation($v, $population);
isa_ok( $ldfcs, 'ARRAY', 'fetch_all_by_Variation returns array ref' );
ok(scalar @$ldfcs == 1, 'fetch_all_by_Variation returns 1 ld feature container');

$ldfc = $ldfcs->[0];
ok($ldfc->isa('Bio::EnsEMBL::Variation::LDFeatureContainer'), "is LDFeatureContainer");

$ld_values = $ldfc->get_all_ld_values;

my $result = {
  rs79944118 => { r2 => 0.087731, d_prime => 0.999965 },
  rs1333047 =>  { r2 => 0.050071, d_prime => 0.999871 },
  rs1333048 =>  { r2 => 0.684916, d_prime => 0.999999 },
  rs10757279 => { r2 => 0.614303, d_prime => 0.817026 },
  rs1333050 =>  { r2 => 0.354847, d_prime => 0.648413 },
  rs4977575 =>  { r2 => 0.050071, d_prime => 0.999871 },
};

foreach my $ld_value (@$ld_values) {
  my $variation_name2 = $ld_value->{variation_name2};
  my $r2 = $ld_value->{r2};
  my $d_prime = $ld_value->{d_prime};
  cmp_ok($result->{$variation_name2}->{r2}, '==', $r2, "fetch_all_by_Variation r2 for $variation_name2 and rs1333049");
  cmp_ok($result->{$variation_name2}->{d_prime}, '==', $d_prime, "fetch_all_by_Variation d_prime for $variation_name2 and rs1333049");
}

#fetch_by_VariationFeature
my $vf = ($v->get_all_VariationFeatures)->[0];
$ldfc = $ldfca->fetch_by_VariationFeature($vf, $population);
foreach my $ld_value (@$ld_values) {
  my $variation_name2 = $ld_value->{variation_name2};
  my $r2 = $ld_value->{r2};
  my $d_prime = $ld_value->{d_prime};
  cmp_ok($result->{$variation_name2}->{r2}, '==', $r2, "fetch_by_VariationFeature r2 for $variation_name2 and rs1333049");
  cmp_ok($result->{$variation_name2}->{d_prime}, '==', $d_prime, "fetch_by_VariationFeature d_prime for $variation_name2 and rs1333049");
}

#fetch_by_VariationFeature, population is not defined
$ldfc = $ldfca->fetch_by_VariationFeature($vf);
foreach my $ld_value (@$ld_values) {
  my $variation_name2 = $ld_value->{variation_name2};
  my $r2 = $ld_value->{r2};
  my $d_prime = $ld_value->{d_prime};
  cmp_ok($result->{$variation_name2}->{r2}, '==', $r2, "fetch_by_VariationFeature population not defined r2 for $variation_name2 and rs1333049");
  cmp_ok($result->{$variation_name2}->{d_prime}, '==', $d_prime, "fetch_by_VariationFeature population not defined d_prime for $variation_name2 and rs1333049");
}

# set max_snp_distance
$ldfca->max_snp_distance(500_000);
cmp_ok($ldfca->max_snp_distance, '==', 500_000, "set/get max_snp_distance");
$ldfc = $ldfca->fetch_by_VariationFeature($vf);
$ld_values = $ldfc->get_all_ld_values;
cmp_ok(scalar @$ld_values, '==', 6, "Number of LD values after changing max_snp_distance");

my $r2 = $ldfca->min_r2;
cmp_ok($r2, '==', 0.0, "set/get min r2");
$r2 = $ldfca->min_r2(0.5);
cmp_ok($r2, '==', 0.5, "set/get min r2");

my $d_prime = $ldfca->min_d_prime;
cmp_ok($d_prime, '==', 0.0, "set/get min d_prime");
$d_prime = $ldfca->min_d_prime(0.9);
cmp_ok($d_prime, '==', 0.9, "set/get min d_prime");

$ldfc = $ldfca->fetch_by_VariationFeature($vf);
$ld_values = $ldfc->get_all_ld_values;

cmp_ok(scalar @$ld_values, '==', 1, "Number of LD values after changing r2 and d_prime values");

# back to default
$r2 = $ldfca->min_r2(0.0);
cmp_ok($r2, '==', 0.0, "set/get min r2");

$d_prime = $ldfca->min_d_prime(0.0);
cmp_ok($d_prime, '==', 0.0, "set/get min d_prime");

#fetch_by_VariationFeatures
my $vf1 = ($va->fetch_by_name('rs1333047')->get_all_VariationFeatures)->[0];
my $vf2 = ($va->fetch_by_name('rs72655407')->get_all_VariationFeatures)->[0];
my $population_id = $population->dbID;
$ldfc = $ldfca->fetch_by_VariationFeatures([$vf1, $vf2], $population);

cmp_ok($ldfc->get_r_square($vf1, $vf2, $population_id), '==', 0.063754, "fetch_by_VariationFeatures with 2 VFs r2");
cmp_ok($ldfc->get_d_prime($vf1, $vf2, $population_id), '==', 0.999996, "fetch_by_VariationFeatures with 2 VFs d_prime");

my $ld_hash = $ldfc->get_all_ld_values(1)->[0];
cmp_ok($ld_hash->{variation_name1}, 'eq', 'rs1333047', "Keep order of input variants in output if running as pairwise calculation");

my $vf3 = ($va->fetch_by_name('rs4977575')->get_all_VariationFeatures)->[0];
$ldfc = $ldfca->fetch_by_VariationFeatures([$vf1, $vf2, $vf3], $population);
cmp_ok($ldfc->get_r_square($vf1, $vf2, $population_id), '==', 0.063754, "fetch_by_VariationFeatures with 3 VFs r2");
cmp_ok($ldfc->get_d_prime($vf1, $vf2, $population_id), '==', 0.999996, "fetch_by_VariationFeatures with 3 VFs d_prime");

cmp_ok($ldfc->get_r_square($vf1, $vf3, $population_id), '==', 1.0, "fetch_by_VariationFeatures with 3 VFs r2");
cmp_ok($ldfc->get_d_prime($vf1, $vf3, $population_id), '==', 1.0, "fetch_by_VariationFeatures with 3 VFs d_prime");

cmp_ok($ldfc->get_r_square($vf2, $vf3, $population_id), '==', 0.063754, "fetch_by_VariationFeatures with 3 VFs r2");
cmp_ok($ldfc->get_d_prime($vf2, $vf3, $population_id), '==', 0.999996, "fetch_by_VariationFeatures with 3 VFs d_prime");

my $populations = $ldfca->get_populations_hash_by_Slice;
cmp_ok(scalar keys %$populations, '==', 48, "Get number of LD populations");

throws_ok { $ldfca->fetch_by_Slice('Slice'); } qr/Bio::EnsEMBL::Slice arg or listref of Bio::EnsEMBL::Slice expected/, 'fetch_all_by_Slice Throw on wrong argument';
throws_ok { $ldfca->fetch_by_Slice(['Slice']); } qr/Bio::EnsEMBL::Slice arg expected/, 'fetch_all_by_Slice Throw on wrong argument';
throws_ok { $ldfca->fetch_by_Slice([$vf1]); } qr/Bio::EnsEMBL::Slice arg expected/, 'fetch_all_by_Slice Throw on wrong argument';

throws_ok { $ldfca->fetch_all_by_Variation('Variation'); } qr/Bio::EnsEMBL::Variation::Variation arg expected/, 'fetch_all_by_Variation Throw on wrong argument';
throws_ok { $ldfca->fetch_by_VariationFeature('VariationFeature'); } qr/Bio::EnsEMBL::Variation::VariationFeature arg expected/, 'fetch_by_VariationFeature Throw on wrong argument';
throws_ok { $ldfca->fetch_by_VariationFeatures(['VariationFeature']); } qr/Bio::EnsEMBL::Variation::VariationFeature arg expected/, 'fetch_by_VariationFeatures Throw on wrong argument';
throws_ok { $ldfca->fetch_by_VariationFeatures('VariationFeature'); } qr/Listref of Bio::EnsEMBL::Variation::VariationFeature args expected/, 'fetch_by_VariationFeatures Throw on wrong argument';

# Test for matching data from VCF file by location and allele instead of variant name

ok($vdb->vcf_config_file($dir.'/ld_no_id_vcf_config.json') eq $dir.'/ld_no_id_vcf_config.json', "DBAdaptor ld_no_id_vcf_config file");
$vca = $vdb->get_VCFCollectionAdaptor();
$coll = $vca->fetch_by_id('ld');
$temp = $coll->filename_template();
$temp =~ s/###t\-root###/$dir/;
$coll->filename_template($temp);
$ldfca->db->use_vcf(1);

$ldfc = $ldfca->fetch_by_Slice($slice, $population);
$ld_values = $ldfc->get_all_ld_values;

foreach my $result (@$results_fetch_by_Slice) {
  my $vf1 = ($va->fetch_by_name($result->{vf1_name})->get_all_VariationFeatures)->[0];
  my $vf2 = ($va->fetch_by_name($result->{vf2_name})->get_all_VariationFeatures)->[0];
  my $population_id = $result->{population_id};
  my $r2 = $result->{r2};
  my $d_prime = $result->{d_prime};
  my $test_name = $result->{test_name};  
  cmp_ok($ldfc->get_r_square($vf1, $vf2, $population_id), '==', $r2, "$test_name r2 no_id_in_vcf");
  cmp_ok($ldfc->get_d_prime($vf1, $vf2, $population_id), '==', $d_prime, "$test_name d_prime no_id_in_vcf");
}

#fetch_all_by_Variation
$v = $va->fetch_by_name('rs1333049');
$ldfcs = $ldfca->fetch_all_by_Variation($v, $population);
isa_ok( $ldfcs, 'ARRAY', 'fetch_all_by_Variation returns array ref' );
ok(scalar @$ldfcs == 1, 'fetch_all_by_Variation returns 1 ld feature container');

$ldfc = $ldfcs->[0];
ok($ldfc->isa('Bio::EnsEMBL::Variation::LDFeatureContainer'), "is LDFeatureContainer");

$ld_values = $ldfc->get_all_ld_values;

foreach my $ld_value (@$ld_values) {
  my $variation_name2 = $ld_value->{variation_name2};
  my $r2 = $ld_value->{r2};
  my $d_prime = $ld_value->{d_prime};
  cmp_ok($result->{$variation_name2}->{r2}, '==', $r2, "fetch_all_by_Variation r2 for $variation_name2 and rs1333049 no_id_in_vcf");
  cmp_ok($result->{$variation_name2}->{d_prime}, '==', $d_prime, "fetch_all_by_Variation d_prime for $variation_name2 and rs1333049 no_id_in_vcf");
}

#fetch_by_VariationFeature
$vf = ($v->get_all_VariationFeatures)->[0];
$ldfc = $ldfca->fetch_by_VariationFeature($vf, $population);
foreach my $ld_value (@$ld_values) {
  my $variation_name2 = $ld_value->{variation_name2};
  my $r2 = $ld_value->{r2};
  my $d_prime = $ld_value->{d_prime};
  cmp_ok($result->{$variation_name2}->{r2}, '==', $r2, "fetch_by_VariationFeature r2 for $variation_name2 and rs1333049 no_id_in_vcf");
  cmp_ok($result->{$variation_name2}->{d_prime}, '==', $d_prime, "fetch_by_VariationFeature d_prime for $variation_name2 and rs1333049 no_id_in_vcf");
}

# test that API can deal with negative start values if max_snp_distance is causing slice start to be negative
$ldfca->max_snp_distance(22_125_504); # use VF start and max_snp_distance to create slice with a negative start
$ldfc = $ldfca->fetch_by_VariationFeature($vf, $population);
cmp_ok(scalar @{$ldfc->get_all_ld_values}, '==', 6, "No error is thrown if max_snp_distance creates negative start value");
# reset max snp distance
$ldfca->max_snp_distance(100_000);

#fetch_by_VariationFeature, population is not defined
$ldfc = $ldfca->fetch_by_VariationFeature($vf);
foreach my $ld_value (@$ld_values) {
  my $variation_name2 = $ld_value->{variation_name2};
  my $r2 = $ld_value->{r2};
  my $d_prime = $ld_value->{d_prime};
  cmp_ok($result->{$variation_name2}->{r2}, '==', $r2, "fetch_by_VariationFeature population not defined r2 for $variation_name2 and rs1333049 no_id_in_vcf");
  cmp_ok($result->{$variation_name2}->{d_prime}, '==', $d_prime, "fetch_by_VariationFeature population not defined d_prime for $variation_name2 and rs1333049 no_id_in_vcf");
}

my $c = Bio::EnsEMBL::Variation::VCFCollection->new(
-id =>  'ld_without_rs_in_vcf_strict_name_match',
-species => 'homo_sapiens',
-assembly => 'GRCh37',
-type =>  'local',
-strict_name_match => 1,
-filename_template => $dir.'/testdata/ld_no_rs_in_vcf.vcf.gz',
-chromosomes => [9],
-sample_prefix => '1000GENOMES:phase_1:',
-adaptor => $vca
);

$vca->remove_VCFCollection_by_ID('ld');
$vca->add_VCFCollection($c);
$ldfc = $ldfca->fetch_by_Slice($slice, $population);
$ld_values = $ldfc->get_all_ld_values;

cmp_ok(scalar @{$ldfc->get_all_ld_values(0)}, '==', 0, "ld_without_rs_in_vcf_strict_name_match match variation feature by name");
cmp_ok(scalar @{$ldfc->get_all_ld_values(1)}, '==', 14, "ld_without_rs_in_vcf_strict_name_match do not match variation feature by name");

$ld_hash = $ldfc->get_all_ld_values(1)->[0];
my $variation1 = $ld_hash->{variation_name1};
my $variation2 = $ld_hash->{variation_name2};
ok($variation1 eq '.' && $variation2 eq '.', 'return name from vcf file');

$c = Bio::EnsEMBL::Variation::VCFCollection->new(
-id =>  'ld_without_rs_in_vcf',
-species => 'homo_sapiens',
-assembly => 'GRCh37',
-type =>  'local',
-filename_template => $dir.'/testdata/ld_no_rs_in_vcf.vcf.gz',
-chromosomes => [9],
-sample_prefix => '1000GENOMES:phase_1:',
-adaptor => $vca
);

$vca->remove_VCFCollection_by_ID('ld_without_rs_in_vcf_strict_name_match');
$vca->add_VCFCollection($c);
delete $ldfca->{_cached};
$ldfc = $ldfca->fetch_by_Slice($slice, $population);
$ld_values = $ldfc->get_all_ld_values;

cmp_ok(scalar @{$ldfc->get_all_ld_values(1)}, '==', 14, "ld_without_rs_in_vcf do not match variation feature by name");
cmp_ok(scalar @{$ldfc->get_all_ld_values(0)}, '==', 14, "ld_without_rs_in_vcf match variation feature by name");

foreach my $result (@$results_fetch_by_Slice) {
  my $vf1 = ($va->fetch_by_name($result->{vf1_name})->get_all_VariationFeatures)->[0];
  my $vf2 = ($va->fetch_by_name($result->{vf2_name})->get_all_VariationFeatures)->[0];
  my $population_id = $result->{population_id};
  my $r2 = $result->{r2};
  my $d_prime = $result->{d_prime};
  my $test_name = $result->{test_name};  
  cmp_ok($ldfc->get_r_square($vf1, $vf2, $population_id), '==', $r2, "$test_name ld_without_rs_in_vcf r2");
  cmp_ok($ldfc->get_d_prime($vf1, $vf2, $population_id), '==', $d_prime, "$test_name ld_without_rs_in_vcf d_prime");
}

$c = Bio::EnsEMBL::Variation::VCFCollection->new(
-id =>  'ld_chr_synonyms',
-species => 'homo_sapiens',
-assembly => 'GRCh37',
-type =>  'local',
-filename_template => $dir.'/testdata/ld_chr_synonyms.vcf.gz',
-use_seq_region_synonyms => 1,
-chromosomes => [9],
-sample_prefix => '1000GENOMES:phase_1:',
-adaptor => $vca
);

$vca->remove_VCFCollection_by_ID('ld_without_rs_in_vcf');
$vca->add_VCFCollection($c);
delete $ldfca->{_cached};
$ldfc = $ldfca->fetch_by_Slice($slice, $population);

cmp_ok(scalar @{$ldfc->get_all_ld_values(1)}, '==', 14, "use chr synonyms -- do not match variation feature by name");
cmp_ok(scalar @{$ldfc->get_all_ld_values(0)}, '==', 14, "use chr synonyms -- match variation feature by name");

done_testing();
