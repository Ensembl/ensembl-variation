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
use Test::More;
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

my $sa = $db->get_SliceAdaptor();
my $vfa = $vdb->get_VariationFeatureAdaptor();
my $va = $vdb->get_VariationAdaptor();
my $pa = $vdb->get_PopulationAdaptor();


my $slice = $sa->fetch_by_region('chromosome', '9', 22124503, 22126503);
#rs4977575 and rs1333050
my $vf1 = $vfa->fetch_by_dbID(3854101);
my $vf2 = $vfa->fetch_by_dbID(1004337);
my $pop_id = 101082;

my $p1 = $pa->fetch_by_name('population');

# fetch_by_Slice
#$ldContainer = $ldfca->fetch_by_Slice($slice, $p1);

my $ld_values;

## VCF
my $dir = $multi->curr_dir();
ok($vdb->vcf_config_file($dir.'/ld_vcf_config.json') eq $dir.'/ld_vcf_config.json', "DBAdaptor vcf_config_file");
my $vca = $vdb->get_VCFCollectionAdaptor();
my $coll = $vca->fetch_by_id('ld');

# now we need to set the filename_template
my $temp = $coll->filename_template();
$temp =~ s/###t\-root###/$dir/;
$coll->filename_template($temp);

# use just VCF
my $p2 = $pa->fetch_by_name('1000GENOMES:phase_1_ASW');
#rs1333047 1004334 rs4977575 3854101
my $vf_a = $vfa->fetch_by_dbID(1004334);
my $vf_b = $vfa->fetch_by_dbID(3854101);

$ldfca->db->use_vcf(2);

$ldContainer = $ldfca->fetch_by_Slice($slice, $p2);
print_container($ldContainer);
$ld_values = count_ld_values($ldContainer);
is($ld_values, 14, "fetch_by_Slice - VCF only");

$ldContainer = $ldfca->fetch_by_VariationFeatures([$vf_a, $vf_b], $p2);
print_container($ldContainer);
$ld_values = count_ld_values($ldContainer);
is($ld_values, 1, "fetch_by_VariationFeatures - VCF only");

# use VCF and DB

$ldfca->db->use_vcf(1);
$ldContainer = $ldfca->fetch_by_Slice($slice);

print_container($ldContainer);
$ld_values = count_ld_values($ldContainer);
is($ld_values, 42, "fetch_by_Slice - VCF and DB");


$ldContainer = $ldfca->fetch_by_VariationFeatures([$vf_a, $vf_b], $p2);
print_container($ldContainer);
$ld_values = count_ld_values($ldContainer);
is($ld_values, 1, "fetch_by_VariationFeatures - VCF and DB");

done_testing();

sub count_ld_values {
  my $container = shift;
  my $ld_values = 0;
  foreach my $key (keys %{$container->{'ldContainer'}}) {
    $ld_values += keys %{$container->{'ldContainer'}->{$key}};
  }
  return $ld_values;
}

sub print_container {
  my $container = shift;
  return if(!$verbose);
  print STDERR "\nContainer name: ", $container->{'name'},"\n";
  foreach my $key (keys %{$container->{'ldContainer'}}) {
    my ($key1,$key2) = split /-/,$key;
    print STDERR "LD values for ", $container->{'variationFeatures'}->{$key1}->variation_name, " and ",$container->{'variationFeatures'}->{$key2}->variation_name;
    foreach my $population (keys %{$container->{'ldContainer'}->{$key}}){
      print STDERR " in population $population:\n d_prime - ",$container->{'ldContainer'}->{$key}->{$population}->{'d_prime'}, "\n r2: ", $container->{'ldContainer'}->{$key}->{$population}->{'r2'}, " \nsample count ",$container->{'ldContainer'}->{$key}->{$population}->{'sample_count'},"\n";
    }
  }
}
