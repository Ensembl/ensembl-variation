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


use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Variation::VariationFeature;

our $verbose = 0;

use_ok('Bio::EnsEMBL::Variation::LDFeatureContainer');
use_ok('Bio::EnsEMBL::Variation::DBSQL::LDFeatureContainerAdaptor');

## examples to be updated & added to test-genome-DBs files


=head
my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();

my $vdb = $multi->get_DBAdaptor('variation');
my $db = $multi->get_DBAdaptor('core');

$vdb->dnadb($db);

my $ldfca = $vdb->get_LDFeatureContainerAdaptor();
my $ldContainer;

ok($ldfca && $ldfca->isa('Bio::EnsEMBL::Variation::DBSQL::LDFeatureContainerAdaptor'));

my $sa = $db->get_SliceAdaptor();

my $slice = $sa->fetch_by_region('chromosome','7');

$ldContainer = $ldfca->fetch_by_Slice($slice,51);

my $ld_values;
print_container($ldContainer);
$ld_values = count_ld_values($ldContainer);
ok($ld_values == 1);

my $vfa = $vdb->get_VariationFeatureAdaptor();

my $vf = $vfa->fetch_by_dbID(153);

$ldContainer = $ldfca->fetch_by_VariationFeature($vf,51);
print_container($ldContainer);
$ld_values = count_ld_values($ldContainer);
ok($ld_values == 1);

sub count_ld_values{
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
=cut

done_testing();
