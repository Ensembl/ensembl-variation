use strict;
use warnings;

use Test::More;
use Bio::EnsEMBL::Variation::Individual;


my $name = 'ind name';
my $description = 'african';
my $gender  = 'Male';
my $mother  = Bio::EnsEMBL::Variation::Individual->new(
    -name => 'mother',
    -description => 'mother',
    -gender  => 'female');

my $father = Bio::EnsEMBL::Variation::Individual->new(
    -name => 'father',
    -description => 'father',
    -gender => 'male');


# test constructor
my $ind = Bio::EnsEMBL::Variation::Individual->new(
    -name => $name,
    -description => $description,
    -gender => $gender,
    -father_individual => $father,
    -mother_individual => $mother);

ok($ind->name() eq $name);
ok($ind->description() eq $description);
ok($ind->gender() eq $gender);
ok($ind->father_Individual() == $father);
ok($ind->mother_Individual() == $mother);
# description
# display
# name
# gender
# type_description
# type_individual
# father_Individual
# mother_Individual
# get_all_child_Individuals
# get_all_Populations

done_testing();

