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
use Bio::EnsEMBL::Variation::VariationSet;
use Bio::EnsEMBL::Variation::DBSQL::VariationSetAdaptor;


my $name          = 'Phenotype-associated variations';
my $description   = 'Variations that have been associated with a phenotype';
my $short_name    = 'ph_variants';


# test constructor
my $variation_set = Bio::EnsEMBL::Variation::VariationSet->new
      ( -dbID        => 12,
        -name        => $name,
        -description => $description,
        -short_name  => $short_name

);



ok($variation_set->name() eq $name, "name");
ok($variation_set->short_name() eq $short_name, "short_name");
ok($variation_set->description() eq $description, "description");



# test getter/setters


ok(test_getter_setter($variation_set, 'name', 'new name'), "get/set new name");
ok(test_getter_setter($variation_set, 'description', 'new description'), "get/set new description");


print "Expecting error message:\n";
# test constructor with a value which is too high

my $variation_set2= Bio::EnsEMBL::Variation::VariationSet->new
      ( -dbID        => 102,
        -name        => $name,
        -description => $description,
        -short_name  => $short_name
      );

done_testing();

