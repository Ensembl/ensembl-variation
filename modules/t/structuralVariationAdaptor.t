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


use Bio::EnsEMBL::Variation::DBSQL::SupportingStructuralVariationAdaptor;
use Bio::EnsEMBL::Variation::DBSQL::VariationSetAdaptor;
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Test::MultiTestDB;

our $verbose = 0;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');

my $vdb = $multi->get_DBAdaptor('variation');
   $vdb->include_failed_variations(1);

my $sva = $vdb->get_StructuralVariationAdaptor();

ok($sva && $sva->isa('Bio::EnsEMBL::Variation::DBSQL::StructuralVariationAdaptor'), "isa sv adaptor");

# test fetch by dbID
print "\n# Test - fetch_by_dbID\n";
my $sv = $sva->fetch_by_dbID(3506221);

ok($sv->variation_name() eq 'esv93078',  'variation name by sv id');
ok($sv->source->name() eq 'DGVa',        'source name by sv id');
ok($sv->study->name() eq 'estd59',       'study name by sv id' );
ok($sv->var_class() eq 'CNV',            'sv class display by sv id' );
ok($sv->class_SO_term() eq 'copy_number_variation', 'sv class SO term by sv id' );
ok($sv->validation_status() eq 'high quality', 'validation status by sv id');

my $clin_sign = $sv->get_all_clinical_significance_states();
ok($clin_sign->[0] eq 'benign' &&
   $clin_sign->[1] eq 'likely benign', 'clinsig by sv id');	

my $study_test = $sv->study;

# test fetch by name
print "\n# Test - fetch_by_name\n";
$sv = $sva->fetch_by_name('esv2421345');

ok($sv->variation_name() eq 'esv2421345', "name by name");
ok($sv->dbID() == 16158654,      "id by name" );
ok($sv->source->name() eq 'DGVa',"source by name");
ok($sv->alias eq 'HM3_CNP_741', "alias by name");

# test store
print "\n# Test - store\n";
delete $sv->{$_} for qw(dbID variation_name);
$sv->variation_name('test');

ok($sva->store($sv), "store");

$sv = $sva->fetch_by_name('test');
ok($sv && $sv->variation_name eq 'test', "fetch stored");


# test fetch by supporting evidence
print "\n# Test - fetch_all_by_supporting_evidence\n";
my $ssva = $vdb->get_SupportingStructuralVariationAdaptor();
my $ssv = $ssva->fetch_by_name('essv194301');
my $sv2 = $sva->fetch_all_by_supporting_evidence($ssv);
ok($sv2->[0]->variation_name() eq 'esv93078', "name by name");

my @sv_dbIDs = (3506221,3506222);
my @sv_names = ('esv93078','esv89107');

# test fetch Iterator by dbID list
print "\n# Test - fetch_Iterator_by_dbID_list\n";
my $sv3 = $sva->fetch_Iterator_by_dbID_list(\@sv_dbIDs);
ok($sv3->next()->variation_name eq $sv_names[0], "iterator by id - 1");
ok($sv3->next()->variation_name eq $sv_names[1], "iterator by id - 2");


## Variation Set ##
my $vsa = $vdb->get_VariationSetAdaptor();
my $vs  = $vsa->fetch_by_name('1000 Genomes - High coverage - Trios');
my @vssv_dbIDs = (3506221,3506222);

# test fetch all dbIDs by VariationSet
print "\n# Test - fetch_all_dbIDs_by_VariationSet\n";
my $sv4 = $sva->fetch_all_dbIDs_by_VariationSet($vs);
ok($sv4->[0] == $vssv_dbIDs[0], "id by VariationSet - 1");
ok($sv4->[1] == $vssv_dbIDs[1], "id by VariationSet - 2");

# test fetch all by VariationSet
print "\n# Test - fetch_all_by_VariationSet\n";
my $sv5 = $sva->fetch_all_by_VariationSet($vs);
ok($sv5->[0]->variation_name eq $sv_names[0], "name by VariationSet - 1");
ok($sv5->[1]->variation_name eq $sv_names[1], "name by VariationSet - 2");

# test fetch Iterator by VariationSet
print "\n# Test - fetch_Iterator_by_VariationSet\n";
my $sv6 = $sva->fetch_Iterator_by_VariationSet($vs);
ok($sv6->next()->variation_name eq $sv_names[0], "iterator by VariationSet - 1");
ok($sv6->next()->variation_name eq $sv_names[1], "iterator by VariationSet - 2");


## Other - BaseStructuralVariationAdaptor ##
my $sv_name = 'CN_674347';

# test fetch all
print "\n# Test - fetch_all\n";
my $sv7= $sva->fetch_all();
ok($sv7->[0]->variation_name eq $sv_name, "sv by all");

# test fetch all somatic
print "\n# Test - fetch_all_somatic\n";
my $sv8= $sva->fetch_all_somatic();
ok($sv8->[0]->variation_name eq 'esv2221103', "sv by all somatic");

# test list dbIDs
print "\n# Test - list_dbIDs\n";
my $sv9= $sva->list_dbIDs();
ok($sv9->[0] == 770260, "sv id by list of dbIDs");

# test fetch by stable_id
print "\n# Test - fetch_by_stable_id\n";
my $sv10= $sva->fetch_by_stable_id($sv_name);
ok($sv10->variation_name eq $sv_name, "sv by stable_id");

# test fetch all by Study
print "\n# Test - fetch_all_by_Study\n";
my $sv11= $sva->fetch_all_by_Study($study_test);
ok($sv11->[0]->variation_name eq $sv_names[0], "sv by study");

# test fetch all by Source
print "\n# Test - fetch_all_by_Source\n";
my $sv12= $sva->fetch_all_by_Source($sv->source);
ok($sv12->[0]->variation_name eq $sv_names[0], "sv by source");

# test get_all_failed_descriptions
print "\n# Test - get_all_failed_descriptions\n";
my $failed_sv = $sva->fetch_by_name('esv902225');
my $sv13= $sva->get_all_failed_descriptions($failed_sv);
ok($sv13->[0] eq 'Variant can not be re-mapped to the current assembly', "failed descriptions");


done_testing();
