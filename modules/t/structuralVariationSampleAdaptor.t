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

our $verbose = 0;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');

my $vdb = $multi->get_DBAdaptor('variation');

my $ssva = $vdb->get_SupportingStructuralVariationAdaptor();
my $svsa = $vdb->get_StructuralVariationSampleAdaptor();

ok($svsa && $svsa->isa('Bio::EnsEMBL::Variation::DBSQL::StructuralVariationSampleAdaptor'), "isa sv sample adaptor");

my $ind_name   = 'NA18635'; 
my $ind_gender = 'Male';
my $ssv_name   = 'essv4067';
my $study_name = 'estd1';


# test fetch by dbID
print "# Test by sv sample id\n";

my $svs = $svsa->fetch_by_dbID(6107305);

ok($svs->individual->name() eq $ind_name,                      'individual name by sv sample id');
ok($svs->individual->gender() eq $ind_gender,                  'individual gender by sv sample id');
ok($svs->structural_variation->variation_name() eq $ssv_name , 'sv name by sv sample id' );
ok($svs->study->name() eq $study_name ,                        'study name by sv sample id' );


# test fetch by ssv object
print "# Test using supporting structural variation object\n";

my $ssv = $ssva->fetch_by_name($ssv_name);
$svs = ($svsa->fetch_all_by_StructuralVariation($ssv))->[0];

ok($svs->individual->name() eq $ind_name,                      'individual name by sv sample id');
ok($svs->individual->gender() eq $ind_gender,                  'individual gender by sv sample id');
ok($svs->structural_variation->variation_name() eq $ssv_name , 'sv name by sv sample id' );
ok($svs->study->name() eq $study_name ,                        'study name by sv sample id' );

done_testing();
