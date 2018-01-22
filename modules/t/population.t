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
use Bio::EnsEMBL::Variation::Population;


my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();

my $vdb = $multi->get_DBAdaptor('variation');

my $pa = $vdb->get_PopulationAdaptor();

# test constructor

my $green_pop = Bio::EnsEMBL::Variation::Population->new
  (-dbID => 123,
   -name => 'Green people',
   -description => 'People who are green',
   -size => 1000,
  );

my $blue_pop = Bio::EnsEMBL::Variation::Population->new
  (-dbID => 124,
   -name => 'Blue people',
   -description => 'People who are blue',
   -size => 1000,
  );


my $sub_pops = [$green_pop, $blue_pop];

my $dbID = 1;
my $name = 'Martians';
my $desc = 'People from the planet Mars';
my $size = 10000;

my $pop = Bio::EnsEMBL::Variation::Population->new
  (-dbID => $dbID,
   -name => $name,
   -description => $desc,
   -size => $size,
   -sub_populations => $sub_pops);

ok($pop->dbID() == $dbID,        'dbID');
ok($pop->name() eq $name,        'name');
ok($pop->description() eq $desc, 'description');
ok($pop->size() == $size,        'size');
ok($pop->get_all_sub_Populations()->[0]->name eq 'Green people', 'sub_pop');

# test getter/setters

ok(test_getter_setter($pop, 'dbID', 123), 'setter/getter dbID');
ok(test_getter_setter($pop, 'name', 'Saturn People'), 'setter/getter name');
ok(test_getter_setter($pop, 'description', 'People from Saturn' ), 'setter/getter description');
ok(test_getter_setter($pop, 'size', 10), 'setter/getter size');



my $purple_pop = Bio::EnsEMBL::Variation::Population->new
  (-dbID => 125,
   -name => 'Purple people',
   -desc => 'People who are purple',
   -size => 1000);

$blue_pop->add_sub_Population($purple_pop);

ok($blue_pop->get_all_sub_Populations()->[0] == $purple_pop, 'add_sub_Population');


# test get_all_super_Populations


my $super_pop = '1000GENOMES:phase_1_EUR';
my $sub_pop   = '1000GENOMES:phase_1_FIN';

# test get all super Populations
my $pop1 = $pa->fetch_by_name($sub_pop);
ok($pop1->get_all_super_Populations->[0]->name eq $super_pop, 'get_all_super_Populations');

# test get all sub Populations
my $pop2 = $pa->fetch_by_name($super_pop);
ok(scalar(@{$pop2->get_all_sub_Populations()}) == 5, 'get_all_sub_Populations');

# test get all synonyms
my $pop3 = $pa->fetch_by_name('PERLEGEN:AFD_EUR_PANEL');
my $synonyms = $pop3->get_all_synonyms();
ok($synonyms->[0] == 1371, 'get_all_synonyms');

# test display group priority
ok($pop2->display_group_priority() == 1, 'display_group_priority');

# test display group name
ok($pop2->display_group_name() eq '1000 Genomes Project Phase 1', 'display_group_name');


done_testing();

