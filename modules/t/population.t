# Copyright 2013 Ensembl
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

BEGIN { $| = 1;
	use Test::More;
	plan tests => 10;
}


use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Variation::Population;

#use Bio::EnsEMBL::Test::MultiTestDB;


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

ok($pop->dbID() == $dbID);
ok($pop->name() eq $name);
ok($pop->description() eq $desc);
ok($pop->size() == $size);
ok($pop->get_all_sub_Populations()->[0]->name eq 'Green people');

# test getter/setters

ok(test_getter_setter($pop, 'dbID', 123));
ok(test_getter_setter($pop, 'name', 'Saturn People'));
ok(test_getter_setter($pop, 'description', 'People from Saturn' ));
ok(test_getter_setter($pop, 'size', 10));



my $purple_pop = Bio::EnsEMBL::Variation::Population->new
  (-dbID => 125,
   -name => 'Purple people',
   -desc => 'People who are purple',
   -size => 1000);

$blue_pop->add_sub_Population($purple_pop);

ok($blue_pop->get_all_sub_Populations()->[0] == $purple_pop);


# test get_all_super_Populations
=head

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();

my $vdb = $multi->get_DBAdaptor('variation');

my $pa = $vdb->get_PopulationAdaptor();

$pop = $pa->fetch_by_name('PACIFIC');

ok(@{$pop->get_all_sub_Populations()} == 25);

$pop = $pop->get_all_sub_Populations->[0];

ok($pop->get_all_super_Populations()->[0]->name() eq 'PACIFIC');

#test get_all_synonyms

$purple_pop->adaptor($pa);
my $synonyms = $purple_pop->get_all_synonyms();

ok($synonyms->[0] == 580);
ok(@{$synonyms} == 1);
=cut
