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
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Variation::Phenotype;


my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();

my $vdb = $multi->get_DBAdaptor('variation');

my $pa = $vdb->get_PhenotypeAdaptor();

# test constructor
my $dbID = 4;
my $stable_id ='testing_stable_id';
my $desc = 'ClinVar: phenotype not specified';
my $class_attrib = 'non_specified';

my $clinvar_pheno = Bio::EnsEMBL::Variation::Phenotype->new
  (-dbID => $dbID,
   -stable_id => $stable_id,
   -description => $desc,
   -class_attrib => $class_attrib);

ok($clinvar_pheno->dbID() == $dbID,                       'dbID');
ok($clinvar_pheno->stable_id() eq $stable_id,             'stable_id');
ok($clinvar_pheno->description() eq $desc,                'description');
ok($clinvar_pheno->class_attrib eq $class_attrib,         'class_attrib');
ok(! defined ($clinvar_pheno->class_attrib_id()),         'class_attrib_id');

# test getter/setters
ok(test_getter_setter($clinvar_pheno, 'dbID', 5), 'setter/getter dbID');
ok(test_getter_setter($clinvar_pheno, 'stable_id', 'Test set stable_id'), 'setter/getter stable_id');
ok(test_getter_setter($clinvar_pheno, 'description', 'Test Setter' ), 'setter/getter description');
ok(test_getter_setter($clinvar_pheno, 'class_attrib_id', 663), 'setter/getter size');

done_testing();
