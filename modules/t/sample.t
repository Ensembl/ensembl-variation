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
use Test::Exception;
use Bio::EnsEMBL::Variation::Individual;
use Bio::EnsEMBL::Variation::Sample;
use Bio::EnsEMBL::Test::MultiTestDB;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');

my $vdb = $multi->get_DBAdaptor('variation');

my $sample_adaptor = $vdb->get_SampleAdaptor;

my $name        = 'ind name';
my $description = 'african';
my $gender      = 'Male';
my $type        = "Outbred";

# test constructor
my $ind = Bio::EnsEMBL::Variation::Individual->new(
  -name => $name,
  -description => $description,
  -gender => $gender,
  -type_individual => $type
);

my $sample_name = 'sample name';
my $sample_description = 'sample desc';
my $display = "DEFAULT";
my $new_display = "DISPLAYABLE";

my $sample = Bio::EnsEMBL::Variation::Sample->new(
  -name => $sample_name,
  -description => $sample_description,
  -display => $display,
  -individual => $ind,
  -adaptor => $sample_adaptor,
);

ok($sample->name() eq $sample_name, "sample name");
ok($sample->description() eq $sample_description, "sample description");
ok($sample->individual->gender() eq $gender, "gender" );
ok($sample->display() eq $display,  "display");
ok($sample->display($new_display) eq $new_display,  "display updated (with argument)"); 
ok($sample->individual->type_individual() eq $type,  "type"); 
ok($sample->has_coverage() == 0,  "default coverage");

$sample->has_coverage(1);
ok($sample->has_coverage() == 1,  "coverage update"); 

throws_ok { $sample->display('WRONG_DISPLAY_VALUE'); } qr/Display flag must be one of/, 'Die on wrong value for display';
throws_ok { $sample->individual('individual'); } qr/Individual argument expected/, 'Die on wrong individual argument';
throws_ok { $sample->study('study'); } qr/Study argument expected/, 'Die on wrong study argument';

$sample->{'study_id'} = 4237;
my $study = $sample->study(); 
ok($study->name eq 'estd1', 'getter for study');
$study = Bio::EnsEMBL::Variation::Study->new(
  -name => 'study_name',
);

$sample->study($study);
$study = $sample->study();
ok($study->name eq 'study_name', 'getter and setter study');

$sample = $sample_adaptor->fetch_by_dbID(101549);
my $populations = $sample->get_all_Populations();
my $concat_populations = join(',', sort map {$_->name} @$populations);
ok($concat_populations eq '1000GENOMES:phase_1_ALL,1000GENOMES:phase_1_CEU,1000GENOMES:phase_1_EUR', 'get_all_Populations');

done_testing();
