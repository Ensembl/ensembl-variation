# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

my $sva = $vdb->get_StructuralVariationAdaptor();

ok($sva && $sva->isa('Bio::EnsEMBL::Variation::DBSQL::StructuralVariationAdaptor'), "isa sv adaptor");

# test fetch by dbID

my $sv = $sva->fetch_by_dbID(3506221);

ok($sv->variation_name() eq 'esv93078',  'variation name by sv id');
ok($sv->source_object->name() eq 'DGVa', 'source name by sv id');
ok($sv->study->name() eq 'estd59',       'study name by sv id' );
ok($sv->var_class() eq 'CNV',            'sv class display by sv id' );
ok($sv->class_SO_term() eq 'copy_number_variation', 'sv class SO term by sv id' );
ok($sv->validation_status() eq 'high quality', 'validation status by sv id');

my $clin_sign = $sv->get_all_clinical_significance_states();
ok($clin_sign->[0] eq 'benign' &&
   $clin_sign->[1] eq 'likely benign', 'clinsig by sv id');	


# test fetch by name
$sv = $sva->fetch_by_name('esv2421345');

ok($sv->variation_name() eq 'esv2421345', "name by name");
ok($sv->dbID() == 16158654,      "id by name" );
ok($sv->source_object->name() eq 'DGVa',"source by name");
ok($sv->alias eq 'HM3_CNP_741', "alias by name");


done_testing();
