# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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
my $display = "DEFAULT";
my $type    = "Outbred";

my $mother  = Bio::EnsEMBL::Variation::Individual->new(
    -name => 'mother',
    -description => 'mother',
    -gender  => 'female');

my $father = Bio::EnsEMBL::Variation::Individual->new(
    -name => 'father',
    -description => 'father',
    -gender => 'male');

my $mother2  = Bio::EnsEMBL::Variation::Individual->new(
    -name => 'mother2',
    -description => 'mother2',
    -gender  => 'female');

my $father2 = Bio::EnsEMBL::Variation::Individual->new(
    -name => 'fathe2r',
    -description => 'father2',
    -gender => 'male');

# test constructor
my $ind = Bio::EnsEMBL::Variation::Individual->new(
    -name => $name,
    -description => $description,
    -gender => $gender,
    -father_individual => $father,
    -mother_individual => $mother,
    -display => $display,
    -type_individual => $type
    );

ok($ind->name() eq $name,   "name");
ok($ind->description() eq $description, " description");
ok($ind->gender() eq $gender, "gender" );
ok($ind->father_Individual() == $father, "father");
ok($ind->mother_Individual() == $mother, "mother" );
ok($ind->display() eq $display,  "display"); 
ok($ind->type_individual() eq $type,  "type"); 
ok($ind->has_coverage() == 0,  "default coverage"); 


my $new_gender  = 'Female';
$ind->gender($new_gender);
ok($ind->gender() eq $new_gender,  "gender update");

$ind->has_coverage(1);
ok($ind->has_coverage() == 1,  "coverage update"); 

$ind->father_Individual($father2);
ok($ind->father_Individual() == $father2, "father update");

$ind->mother_Individual($mother2);
ok($ind->mother_Individual() == $mother2, "mother update " );


# get_all_child_Individuals
# get_all_Populations



done_testing();

