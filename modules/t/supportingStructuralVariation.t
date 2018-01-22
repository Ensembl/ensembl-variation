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
use Test::More;
use Bio::EnsEMBL::Test::MultiTestDB;

use FindBin qw($Bin);

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Variation::SupportingStructuralVariation;
use Bio::EnsEMBL::Variation::Study;
use Bio::EnsEMBL::Variation::Source;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');

my $vdb = $multi->get_DBAdaptor('variation');

my $ssv_adaptor = $vdb->get_SupportingStructuralVariationAdaptor;


# test constructor

### need Source object 
my $source_name        = 'DGVa';
my $source_version     = 201310;
my $source_description = 'Database of Genomic Variants Archive';
my $source_id          = 11;

my $source = Bio::EnsEMBL::Variation::Source->new
  (-dbID        => $source_id,
   -name        => $source_name,
   -version     => $source_version,
   -description => $source_description
);

## need Study object  
my $study_name = 'estd1';
my $study_url = 'ftp://ftp.ebi.ac.uk/pub/databases/dgva/estd1_Redon_et_al_2006';
my $study_description = 'Redon 2006 "Global variation in copy number in the human genome." PMID:17122850 [remapped from build NCBI35]';

my $study = Bio::EnsEMBL::Variation::Study->new
  (-name         => $study_name,
   -url          => $study_url,
   -description  => $study_description,
   -_source_id   => $source->dbID
);


## need Supporting Structural Variantion object 
my $ssv_id    = 18294456;
my $ssv_name  = 'essv4067';
my $ssv_alias = 'NA18635_Chr16_2';

my $ssv = Bio::EnsEMBL::Variation::SupportingStructuralVariation->new
  (-dbID           => $ssv_id,
   -variation_name => $ssv_name,
   -alias          => $ssv_alias,
   -study          => $study,
   -source         => $source,
   -is_evidence    => 1,
);
  
ok($ssv->dbID() eq $ssv_id,             "dbID");
ok($ssv->variation_name() eq $ssv_name, "name");
ok($ssv->display_id() eq $ssv_name,     "display");
ok($ssv->alias() eq $ssv_alias,         "alias");
ok($ssv->study->name() eq $study_name,  "study");
ok($ssv->is_evidence() == 1,            "is_evidence");
ok($ssv->is_somatic() == 0,             "is_somatic");
# source
ok($ssv->source->name() eq $source_name,            'ssv -> source' );
ok($ssv->source_name eq $source_name,               'ssv -> source_name');
ok($ssv->source_description eq $source_description, 'ssv -> source_description');
ok($ssv->source_version eq $source_version,         'ssv -> source_version');


# test getter/setters

ok(test_getter_setter($ssv, 'variation_name', 'newname'), "get/set name");
ok(test_getter_setter($ssv, 'alias','new_alias_1'), "get/set alias");
ok(test_getter_setter($ssv, 'validation_status','new satus'), "get/set validation status");


# test get_all_StructuralVariations

my $real_ssv = $ssv_adaptor->fetch_by_dbID($ssv_id);

my @svs = sort {$a->dbID() <=> $b->dbID()} @{$real_ssv->get_all_StructuralVariations()};

ok($svs[0]->variation_name() eq 'esv2758415',           "sv name");
ok($svs[0]->class_SO_term() eq 'copy_number_variation', "sv SO term");


done_testing();
