use strict;
use warnings;

BEGIN { $| = 1;
	use Test;
	plan tests => 15;
}


use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Variation::Population;

use Bio::EnsEMBL::Test::MultiTestDB;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();

my $vdb = $multi->get_DBAdaptor('variation');


# test constructor

my $green_pop = Bio::EnsEMBL::Variation::Population->new
  (-dbID => 123,
   -name => 'Green people',
   -description => 'People who are green',
   -size => 1000,
   -is_strain => 0);

my $blue_pop = Bio::EnsEMBL::Variation::Population->new
  (-dbID => 124,
   -name => 'Blue people',
   -description => 'People who are blue',
   -size => 1000,
   -is_strain => 1);


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
   -is_strain => 0,
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
ok(test_getter_setter($pop,'is_strain',0));



my $purple_pop = Bio::EnsEMBL::Variation::Population->new
  (-dbID => 125,
   -name => 'Purple people',
   -desc => 'People who are purple',
   -size => 1000);

$blue_pop->add_sub_Population($purple_pop);

ok($blue_pop->get_all_sub_Populations()->[0] == $purple_pop);


# test get_all_super_Populations

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
