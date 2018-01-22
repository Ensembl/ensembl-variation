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


my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');

my $vdb = $multi->get_DBAdaptor('variation');

my $vs_adaptor = $vdb->get_VariationSetAdaptor();
my $vs = $vs_adaptor->fetch_by_short_name('1kg_com');

ok($vs && $vs->isa('Bio::EnsEMBL::Variation::VariationSet'), "isa var set");

ok($vs->short_name() eq '1kg_com',  "variation set short name");
ok($vs->name() eq '1000 Genomes - All - common',  "variation set name");
ok($vs->description() eq 'Variants genotyped by the 1000 Genomes project (phase 1) with frequency of at least 1%',  "variation set description");


my $subsets = $vs->get_all_sub_VariationSets();
ok(scalar(@{$subsets}) == 4,  "count subsets");

my $supersets = $subsets->[0]->get_all_super_VariationSets(); 
ok(scalar(@{$supersets}) == 1,  "count supersets");
ok( $supersets->[0]->short_name() eq '1kg_com',  "super set short name");


## check iterator
my $limit = 10;
my $fetched = 0;
my $it = $vs->get_Variation_Iterator();

ok($it && $it->isa('Bio::EnsEMBL::Utils::Iterator'), "isa var iterator");

my $set_var;
while ($fetched < $limit && $it->has_next()) {
    
    $set_var = $it->next();
    $fetched++;
} 
ok($fetched ==10 ,    "variation iterator limit");
ok($set_var->name() eq 'rs142276873', "last set member is as expected");


## different results if fails are included
$vs_adaptor->db->include_failed_variations(1);
my $vsf = $vs_adaptor->fetch_by_short_name('1kg_com');
my $itf = $vsf->get_Variation_Iterator();

my $fset_var;
my $f_fetched = 0;
while ($f_fetched < $limit && $itf->has_next()) {
    $fset_var = $itf->next();
    $f_fetched++;
} 
ok($f_fetched ==10 ,    "variation (inc fail) iterator limit");
ok($fset_var->name() eq 'rs117161559', "last (inc fail) set member is as expected");


done_testing();
