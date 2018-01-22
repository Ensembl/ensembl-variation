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
use Data::Dumper;
use Test::Exception;
use Test::More;
use Bio::EnsEMBL::Test::MultiTestDB;

use FindBin qw($Bin);

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Variation::Source;
use Bio::EnsEMBL::Variation::Study;
use Bio::EnsEMBL::Variation::Sample;
use Bio::EnsEMBL::Variation::Individual;
use Bio::EnsEMBL::Variation::SupportingStructuralVariation;
use Bio::EnsEMBL::Variation::StructuralVariationSample;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');

my $vdb  = $multi->get_DBAdaptor('variation');

my $svs_adaptor = $vdb->get_StructuralVariationSampleAdaptor;


# test constructor

## need Source object 
my $source_name = 'DGVa';
my $source_id   = 11;

my $source = Bio::EnsEMBL::Variation::Source->new
  (-dbID => $source_id,
   -name => $source_name
);

## need Study object 
my $study_name = 'estd1';

my $study = Bio::EnsEMBL::Variation::Study->new
  (-name       => $study_name,
   -_source_id => $source->dbID
);

## need strain/individual (Individual object)

my $ind_id     = 2105;
my $ind_name   = 'NA18635';
my $ind_gender = 'Male';

my $ind = Bio::EnsEMBL::Variation::Individual->new
  (-dbID        => $ind_id,
   -name        => $ind_name,
   -gender      => $ind_gender
);

## need Sample object 

my $sample_id     = 8675;
my $sample_name   = 'NA19122';

my $sample = Bio::EnsEMBL::Variation::Sample->new
  (-dbID => $sample_id,
   -name => $sample_name,
   -individual => $ind
);


## need Supporting Structural Variantion object 
my $ssv_id   = 18294456;
my $ssv_name = 'essv4067';

my $ssv = Bio::EnsEMBL::Variation::SupportingStructuralVariation->new
  (-dbID   => $ssv_id,
   -name   => $ssv_name,
   -study  => $study,
   -source => $source
);

my $dbID     = 6107305;
my $zygosity = 1;

my $svs = Bio::EnsEMBL::Variation::StructuralVariationSample->new
  (-dbID                     => $dbID,
   -_structural_variation_id => $ssv->dbID,
   -_strain_id               => $ind->dbID,
   -strain                   => $ind,
   -sample                   => $sample,
   -study                    => $study,
   -zygosity                 => $zygosity,
   -adaptor                  => $svs_adaptor
  );

ok($svs->dbID() eq $dbID,                                     'dbID');
ok($svs->structural_variation->variation_name() eq $ssv_name, 'ssv name');
# Strain
ok($svs->strain->name() eq $ind_name,                         'strain name');
ok($svs->strain->gender() eq $ind_gender,                     'strain gender ');
# Individual
ok($svs->sample->individual->name() eq $ind_name,             'individual name');
ok($svs->sample->individual->gender() eq $ind_gender,         'individual gender ');
# Sample
ok($svs->sample->name() eq $sample_name,                      'sample name');
# Study
ok($svs->study->name() eq $study_name ,                       'study name' );
# Zygosity
ok($svs->zygosity() == $zygosity ,                            'zygosity' );


# test structural variation object
my $sv = $svs->structural_variation();
ok($svs->structural_variation($sv), 'structural_variation object (using argument)');

# test study object
my $svs_study = $svs->study();
ok($svs->study($svs_study), 'study object (using argument)');

my $svs_hash = {
  dbID                     => $dbID,
  _structural_variation_id => $ssv->dbID,
  _strain_id               => $ind->dbID,
  sample                   => $sample,
  study                    => $study,
  zygosity                 => $zygosity,
  adaptor                  => $svs_adaptor
};
$svs = Bio::EnsEMBL::Variation::StructuralVariationSample->new_fast($svs_hash);
ok($svs->dbID() eq $dbID, 'dbID');

throws_ok { $svs->structural_variation('structural_variation'); } qr/SupportingStructuralVariation argument expected/, 'Throw on wrong argument structural_variation';
throws_ok { $svs->study('study'); } qr/Study argument expected/, 'Throw on wrong argument study';

my $individual = $svs->strain();
ok($individual->name eq 'NA18635', 'get individual');

throws_ok { $svs->strain('strain'); } qr/Individual argument expected/, 'Throw on wrong argument strain';

my $strain = $svs->strain();
ok($strain->name() eq 'NA18635', 'get strain');

done_testing();
