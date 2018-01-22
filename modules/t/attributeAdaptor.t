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
use Data::Dumper;

use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Test::MultiTestDB;


our $verbose = 0;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');

my $vdb = $multi->get_DBAdaptor('variation');

my $aa = $vdb->get_AttributeAdaptor();

ok($aa && $aa->isa('Bio::EnsEMBL::Variation::DBSQL::AttributeAdaptor'), "Is an attribute adaptor");


# Attrib type
my $attrib_type_id = 1;
my $attrib_type_code = 'SO_accession';
ok($aa->attrib_id_for_type_code($attrib_type_code) eq $attrib_type_id, 'attrib_type_id by attrib type "code"');
ok($aa->attrib_type_code_for_attrib_type_id($attrib_type_id) eq $attrib_type_code, 'Attrib type "code" by attrib_type_id');

my $attrib_type_code2 = 'short_name';
my $attrib_type_name  = 'Short name';
ok($aa->attrib_type_name_for_attrib_type_code($attrib_type_code2) eq $attrib_type_name, 'Attrib type "name" by attrib type "code"');


# Attrib
my $attrib_id    = 3; 
my $attrib_value = 'SNP';
my $attrib_type_code3 = 'display_term';
ok($aa->attrib_value_for_id($attrib_id) eq $attrib_value,                         'Attrib "value" by attrib_id');
ok($aa->attrib_id_for_type_value($attrib_type_code3,$attrib_value) eq $attrib_id, 'attrib_id by attrib type "code" and attrib "value"'); 

# Sequence ontology
my $SO_term = 'missense_variant';
my $SO_accession = 'SO:0001583';
my $display_term = 'NON_SYNONYMOUS_CODING';
ok($aa->display_term_for_SO_term($SO_term) eq $display_term,           'Display term by SO term');
ok($aa->SO_accession_for_SO_term($SO_term) eq $SO_accession,           'SO accession by SO term');
ok($aa->SO_term_for_SO_accession($SO_accession) eq $SO_term,           'SO term by SO accession');
ok($aa->display_term_for_SO_accession($SO_accession) eq $display_term, 'Display term by SO accession');

done_testing();
