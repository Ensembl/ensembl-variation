use strict;
use warnings;

BEGIN { $| = 1;
	use Test;
	plan tests => 10;
}


use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Variation::Individual;


my $name = 'ind name';
my $description = 'african';
my $gender  = 'Male';
my $mother  = Bio::EnsEMBL::Variation::Individual->new
  (-name => 'mother',
   -description => 'mother',
   -gender  => 'female');

my $father = Bio::EnsEMBL::Variation::Individual->new
  (-name => 'father',
   -description => 'father',
   -gender => 'male');


# test constructor
my $ind = Bio::EnsEMBL::Variation::Individual->new
  (-name => $name,
   -description => $description,
   -gender => $gender,
   -father_individual => $father,
   -mother_individual => $mother);


ok($ind->name() eq $name);
ok($ind->description() eq $description);
ok($ind->gender() eq $gender);
ok($ind->father_Individual() == $father);
ok($ind->mother_Individual() == $mother);


#
# test getter/setters
#

ok(test_getter_setter($ind, 'name', 'new name'));
ok(test_getter_setter($ind, 'description', 'new description'));
ok(test_getter_setter($ind, 'gender', 'Female'));


$mother  = Bio::EnsEMBL::Variation::Individual->new
  (-name => 'new mother',
   -description => 'new mother',
   -gender  => 'female');

$father = Bio::EnsEMBL::Variation::Individual->new
  (-name => 'new father',
   -description => 'new father',
   -gender => 'male');

ok(test_getter_setter($ind, 'mother_Individual', $mother));
ok(test_getter_setter($ind, 'father_Individual', $father));



### TODO: add test for get_all_child_Individuals once implemented
