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

