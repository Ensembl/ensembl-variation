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

use lib 't';

use strict;
use warnings;

use Test::More;

use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Test::MultiTestDB;

use Data::Dumper;
use Bio::EnsEMBL::Variation::LDFeatureContainer;
use Bio::EnsEMBL::Variation::VariationFeature;
use Bio::EnsEMBL::Variation::Variation;

our $verbose = 0;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');
my $vdb = $multi->get_DBAdaptor('variation');
my $ldfca = $vdb->get_LDFeatureContainerAdaptor;

#test constructor
my $v1 = Bio::EnsEMBL::Variation::Variation->new(-name => 'rs193', -source => 'dbSNP');
my $v2 = Bio::EnsEMBL::Variation::Variation->new(-name => 'rs203', -source => 'dbSNP');
my $v3 = Bio::EnsEMBL::Variation::Variation->new(-name => 'rs848', -source => 'dbSNP');
my $v4 = Bio::EnsEMBL::Variation::Variation->new(-name => 'rs847', -source => 'dbSNP');


my $vf1 = Bio::EnsEMBL::Variation::VariationFeature->new(
	-dbID => 153,
	-start => 27686081,
	-end => 27686081,
	-strand => 1,
	-variation_name => 'rs193',
	-map_weight => 1,
	-allele_string => 'C/T',
	-variation => $v1
);

my $vf2 = Bio::EnsEMBL::Variation::VariationFeature->new(
	-dbID => 163,
	-start => 27689871,
	-end => 27689871,
	-strand => 1,
	-variation_name => 'rs203',
	-map_weight => 1,
	-allele_string => 'T/C',
	-variation => $v2
);

my $vf3 = Bio::EnsEMBL::Variation::VariationFeature->new(
	-dbID => 749,
	-start => 132072716,
	-end => 132072716,
	-strand => -1,
	-variation_name => 'rs848',
	-map_weight => 1,
	-allele_string => 'T/G',
	-variation => $v3
);

my $vf4 = Bio::EnsEMBL::Variation::VariationFeature->new(
	-dbID => 748,
	-start => 132072885,
	-end => 132072885,
	-strand => -1,
	-variation_name => 'rs847',
	-map_weight => 1,
	-allele_string => 'A/G',
	-variation => $v4
);

my $ldContainer = Bio::EnsEMBL::Variation::LDFeatureContainer->new(
  '-adaptor' => $ldfca,
	'-name' => 'container_1',
	'-ldContainer' => { 
		'27686081-27689871' => { 
			51 => {
				'd_prime'            => 0.533013,
				'r2'                 => 0.258275,
				'sample_count'       => 42
			},
			140 => {
				'd_prime'            => 0.999887,
				'r2'                 => 0.642712,
				'sample_count'       => 10
			}
		},
		'132072716-132072885' => {
			140 => {
				'd_prime'            => 0.999924,
				'r2'                 => 0.312452,
				'sample_count'       => 22
			}
		}
	},
	'-pos2name' => {
		27686081  => 'rs193',
		27689871  => 'rs203',
		132072716 => 'rs848',
		132072885 => 'rs847'
	},
	'-pos2vf' => {
		27686081  => $vf1,
		27689871  => $vf2,
		132072716 => $vf3,
		132072885 => $vf4
	}
);

ok($ldContainer->name()  eq 'container_1', "container name");

print_container($ldContainer);

# test getter_setter

ok(test_getter_setter($ldContainer,'name','container_new_name'), "name setter" );


#test methods
my $variations = $ldContainer->get_variations();
ok(@{$variations} == 4, "count variants in container");

#to check how to get the r_square value for 2 variation_features with a known and an unknown population
my $r_square;
$r_square = $ldContainer->get_r_square($vf1,$vf2,51);
ok($r_square == 0.258275, "r squared for vfs and population");

$r_square = $ldContainer->get_r_square($vf1,$vf2);
ok($r_square == 0.642712, "r squared for vfs, no population");


#to check how to get the d_prime value for 2 variation_features with a known and an unknown population
my $d_prime;
$d_prime = $ldContainer->get_d_prime($vf3,$vf4,140);
ok($d_prime == 0.999924, "d prime for vfs, with population");

$d_prime = $ldContainer->get_d_prime($vf1,$vf2);
ok($d_prime == 0.999887, "d prime for vfs, no population");

#check method to get ALL ld values in container (d_prime, r2, and sample_count
my $ld_values;
$ld_values = $ldContainer->get_all_ld_values();
ok(@{$ld_values} == 2, "count total stats");
my $r_squares = $ldContainer->get_all_r_square_values();
ok(@{$r_squares} == 2, "count r squared");
my $d_primes = $ldContainer->get_all_d_prime_values();
ok(@{$d_primes} == 2, "count d prime");

#check method to retrieve populations in a container
my $populations = $ldContainer->get_all_populations();
ok(@{$populations} == 2, "count populations");
$populations = $ldContainer->get_all_populations($vf3,$vf4);
ok($populations->[0] == 140, "population id");

done_testing();

sub print_container {
  my $container = shift;
  return if(!$verbose);
 
  print STDERR "\nContainer name: ", $container->{'name'},"\n";
  foreach my $key (keys %{$container->{'ldContainer'}}) {
      my ($key1,$key2) = split /-/,$key;
      print STDERR "LD values for ", $container->{'variationFeatures'}->{$key1}->variation_name, " and ",$container->{'variationFeatures'}->{$key2}->variation_name;
      foreach my $population (keys %{$container->{'ldContainer'}->{$key}}){
	  print STDERR " in population $population:\n d_prime - ",$container->{'ldContainer'}->{$key}->{$population}->{'d_prime'}, "\n r2: ", $container->{'ldContainer'}->{$key}->{$population}->{'r2'}, "\n snp_distance: ",$container->{'ldContainer'}->{$key}->{$population}->{'snp_distance_count'}, " \nsample count ",$container->{'ldContainer'}->{$key}->{$population}->{'sample_count'},"\n";
      }
  }

}
