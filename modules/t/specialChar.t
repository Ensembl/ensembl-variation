#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2025] EMBL-European Bioinformatics Institute
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
use FindBin qw($Bin);

use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Variation::Utils::SpecialChar qw(replace_char replace_hex decode_text);


my $oldString = "Hôpitaux de Paris";
my $newString = replace_char($oldString);
ok($newString eq "Hopitaux de Paris", "special char updated");

$oldString = "Auto-Immunite Rhumatoïde,Hôpitaux";
$newString = replace_char($oldString);
ok($newString eq "Auto-Immunite Rhumatoide,Hopitaux", "special char(ï, ô) updated");

$oldString = "Sjögren-Larsson syndrome";
$newString = replace_char($oldString);
ok($newString eq "Sjogren-Larsson syndrome", "special char(ö) updated");

$oldString = "bmi < 30";
$newString = replace_char($oldString);
ok($newString eq "bmi less than 30", "special char(<) updated");

$oldString = "bmi > 30";
$newString = replace_char($oldString);
ok($newString eq "bmi more than 30", "special char(>) updated");

$oldString = '3-@METHYLGLUTACONIC ACIDURIA, TYPE I';
$newString = replace_char($oldString);
ok($newString eq "3-METHYLGLUTACONIC ACIDURIA, TYPE I", "special char(@) removed");

$oldString = "isom&#xe8;res";
$newString = replace_hex($oldString);
ok($newString eq "isomeres", "hex replaced");

$oldString = "Marinesco-Sj+¦gren syndrome";
$newString = decode_text($oldString);
ok($newString eq "Marinesco-Sjogren syndrome", "text replaced");

done_testing();
