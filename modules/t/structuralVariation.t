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
use Bio::EnsEMBL::Variation::StructuralVariation;
use Bio::EnsEMBL::Variation::Study;
use Bio::EnsEMBL::Variation::Source;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');

my $vdb  = $multi->get_DBAdaptor('variation');

my $sv_adaptor = $vdb->get_StructuralVariationAdaptor;


# test constructor

## need Source object 
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
my $study_name = 'estd59';
my $study_url = 'ftp://ftp.ebi.ac.uk/pub/databases/dgva/estd59_1000_Genomes_Consortium_Pilot_Project';
my $study_description = '1000 Genomes Project Consortium - Pilot Project. PMID:20981092';

my $study = Bio::EnsEMBL::Variation::Study->new
  (-name         => $study_name,
   -url          => $study_url,
   -description  => $study_description,
   -_source_id   => $source->dbID
);


my $dbID = 3506221;
my $name = 'esv93078';
my $alias = 'alias_1';
my $class_SO_term = 'copy_number_variation';
my $SO_accession  = 'SO:0001019';
my $validation_status = 'high quality';
my $clin_signs = ['benign','likely benign'];


my $sv = Bio::EnsEMBL::Variation::StructuralVariation->new
  (-dbID                  => $dbID,
   -variation_name        => $name,
   -source                => $source,
   -study                 => $study,
   -class_SO_term         => $class_SO_term,
   -alias                 => $alias,
   -adaptor               => $sv_adaptor,
   -validation_status     => $validation_status,
   -clinical_significance => $clin_signs,
   -is_evidence           => 0,
   -is_somatic            => 0,
  );

ok($sv->dbID() eq $dbID,                           "dbID");
ok($sv->variation_name() eq $name,                 "name");
ok($sv->display_id() eq $name,                     "display");
ok($sv->alias() eq $alias,                         "alias");
ok($sv->study->name() eq $study_name,              "study");
ok($sv->class_SO_term() eq $class_SO_term,         "class SO term" );
ok($sv->class_SO_accession() eq $SO_accession,     "class SO accesssion" );
ok($sv->validation_status() eq $validation_status, "validation status");
ok($sv->is_evidence() eq 0,                        "is_evidence");
ok($sv->is_somatic() eq 0,                         "is_somatic");
# source
ok($sv->source->name() eq $source_name,     'sv -> source' );
ok($sv->source_name eq $source_name,               'sv -> source_name');
ok($sv->source_description eq $source_description, 'sv -> source_description');
ok($sv->source_version eq $source_version,         'sv -> source_version');



my $clin_sign = $sv->get_all_clinical_significance_states();
ok($clin_sign->[0] eq 'benign' &&
   $clin_sign->[1] eq 'likely benign', 'clinsig by sv id');	


# test getter/setters

ok(test_getter_setter($sv, 'variation_name', 'newname'), "get/set name");
ok(test_getter_setter($sv, 'alias','new_alias_1'), "get/set alias");
ok(test_getter_setter($sv, 'validation_status','new satus'), "get/set validation status");


#test get_all_SupportingStructuralVariants

my $real_sv = $sv_adaptor->fetch_by_dbID($dbID);

my @ssvs = sort {$a->dbID() <=> $b->dbID()} @{$real_sv->get_all_SupportingStructuralVariants()};

ok(@ssvs == 2,                                      "ssv count" );
ok($ssvs[0]->variation_name() eq 'essv194300',      "ssv name");
ok($ssvs[0]->class_SO_term() eq 'copy_number_loss', "ssv SO term");
my $ssv_cs = $ssvs[0]->get_all_clinical_significance_states();
ok($ssv_cs->[0] eq $clin_sign->[0],                 "ssv clin sign");
ok($ssvs[0]->alias eq 'alias_2',                    "ssv alias");



my $sv2 = $sv_adaptor->fetch_by_name($name);


# test source object
my $sv_source = $sv2->source();
ok($sv2->source($sv_source), 'source (using argument)');

# test study object
my $sv_study = $sv2->study();
ok($sv2->study($sv_study), 'study object (using argument)');

# test get all StructuralVariationFeatures
my $svfs = $sv2->get_all_StructuralVariationFeatures();
ok($svfs->[0]->seq_region_name eq '8' && $svfs->[0]->seq_region_start == 7823440, 'sv -> get_all_StructuralVariationFeatures');

# test get all PhenotypeFeatures
my $pfs = $sv2->get_all_PhenotypeFeatures();
ok($pfs->[0]->phenotype->description eq 'ACHONDROPLASIA', 'sv -> get_all_PhenotypeFeatures');

# test get all StructuralVariationSamples
my $svss = $sv2->get_all_StructuralVariationSamples();
ok($svss->[0]->sample->name eq 'NA12891', 'sv -> get_all_StructuralVariationSamples');


ok($sv2->summary_as_hash->{'display_id'} eq $name, 'sv-> summary');


done_testing();
