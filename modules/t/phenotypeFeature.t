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
use Data::Dumper;
use Bio::EnsEMBL::Test::MultiTestDB;

use Bio::EnsEMBL::Variation::Source;
use Bio::EnsEMBL::Variation::Study;
use Bio::EnsEMBL::Variation::Variation;
use Bio::EnsEMBL::Variation::Phenotype;
use Bio::EnsEMBL::Variation::PhenotypeFeature;
use Bio::EnsEMBL::Slice;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');

my $vdb = $multi->get_DBAdaptor('variation');
my $db  = $multi->get_DBAdaptor('core');

my $pf_adaptor   = $vdb->get_PhenotypeFeatureAdaptor;


## need source object 
my $source_id             = 1;
my $source_name           = 'dbSNP';
my $source_version        = 138;
my $source_description    = 'Variants (including SNPs and indels) imported from dbSNP';

my $source = Bio::EnsEMBL::Variation::Source->new
  (-dbID           => $source_id,
   -name           => $source_name,
   -version        => $source_version,
   -description    => $source_description
);

## need a study
my $study_name        = "Beverage impact"; 
my $study_description = "Beverage impact consortium";
my $study_url         = 'http://bic.org';
my $study_xref        = 'pubmed/10000'; 
my $study = Bio::EnsEMBL::Variation::Study->new
  (-name               => $study_name,
   -description        => $study_description,
   -source             => $source,
   -url                => $study_url,
   -external_reference => $study_xref
);


my $v_name = 'rs142276873';

## need a feature
my $variation = Bio::EnsEMBL::Variation::Variation->new(-name   => $v_name,
                                                        -source => $source);


## need a slice
my $sa = $db->get_SliceAdaptor();
my $slice = $sa->fetch_by_region('chromosome', '18');


my $external_id  = 12345;
my $p_value      = 0.0000023;
my $risk_allele  = 'G';
my $phenotype_id = 1;
my $desc         = 'Tea Consumption';
my $gene         = 'TEA1';
my $clinsig      = 'protective';
my $or           = 6;
my $beta         = 2;
my $allele_symbol = 't_1';
my $allele_accession = 't_1.1';
my $var_name = "rs12345";
my $external_ref = "RV123";
my @ontology_accessions = ('EFO:00001');

my $phenotype = Bio::EnsEMBL::Variation::Phenotype->new(-DESCRIPTION => $desc, -dbID => $phenotype_id);

my $pf = Bio::EnsEMBL::Variation::PhenotypeFeature->new(
    -slice     => $slice,
    -start     => 23821095,
    -end       => 23821095,
    -phenotype => $phenotype,
    -type      => 'Variation',
    -object    => $variation,
    -source    => $source,
    -study     => $study,
    -is_significant  => 1,
    -attribs   => {
      p_value         => $p_value,
      beta_coef       => $beta,
      odds_ratio      => $or,
      external_id     => $external_id,
      risk_allele     => $risk_allele,
      associated_gene => $gene, 
      clinvar_clin_sig     => $clinsig,
      allele_symbol        => $allele_symbol,     
      allele_accession_id  => $allele_accession,
      variation_names      => $var_name,
      external_reference   => $external_ref 
    },
    -ontology_accessions => \@ontology_accessions,
    );


ok($pf->start() == 23821095,                       "start");
ok($pf->end()   == 23821095,                       "end") ;
ok($pf->external_id() eq $external_id,             "external_id");
ok($pf->risk_allele() eq $risk_allele,             "risk_allele");
ok($pf->p_value ()    eq  $p_value,                "p_value");
ok($pf->beta_coefficient() eq $beta,               "beta_coefficient");
ok($pf->odds_ratio()    eq $or,                    "or");
ok($pf->is_significant()  eq  1,                   "is significant");
ok($pf->clinical_significance() eq $clinsig,       "clinical significance");
ok($pf->type() eq 'Variation',                     "type");
ok($pf->associated_gene  eq $gene,                 "associated_gene");
ok($pf->allele_symbol() eq $allele_symbol,         "allele_symbol");
ok($pf->allele_accession_id() eq $allele_accession,"allele_accession");
ok($pf->source_name() eq $source_name,             "source name");
ok($pf->source_version() eq $source_version,       "source version");
ok($pf->study_name() eq $study_name,               "study name");
ok($pf->study_url()  eq $study_url,                "study_url");
ok($pf->study_description() eq $study_description, "study description");
ok($pf->variation_names() eq $var_name,            "assoc var names");
ok($pf->phenotype() eq $phenotype,                 "phenotype object");
ok($pf->phenotype_id() eq $phenotype_id,           "phenotype ID");
ok($pf->phenotype()->description eq $desc,         "phenotype");
ok($pf->object()->name() eq $v_name,               "variation name");
ok($pf->external_reference() eq $study_xref,       "external reference");
ok($pf->get_all_ontology_accessions->[0] eq $ontology_accessions[0], "ontology accessions");

# test source object
my $pf_object = $pf->object();
ok($pf->object($pf_object)->name eq $v_name, 'pf object (using argument)');

##set var
my $variation2 = Bio::EnsEMBL::Variation::Variation->new(-name   => 'rs1234',
                                                        -source => $source);
$pf->variation($variation2);
ok($pf->variation()->name()   eq 'rs1234',   "updated variation name");

## new source object 
my $source_name2           = 'ClinVar';
my $source_version2        = 138;


my $source2 = Bio::EnsEMBL::Variation::Source->new
  (-name           => $source_name2,
   -version        => $source_version2,
);
$pf->source($source2);
ok($pf->source_name() eq $source_name2,             "update source");


# Tests to populate empty variables ('object', 'variation', 'source')
my $pf2 = Bio::EnsEMBL::Variation::PhenotypeFeature->new(
    -adaptor    => $pf_adaptor,
    -type       => 'Variation',
    -_object_id => $variation->name,
    -_source_id => $source_id,
    -attribs    => {strain_id  => 1592}, # Individual 'NA12891'
   );
ok($pf2->object()->name eq $v_name,    "pf 'object' 2 (using argument)");
ok($pf2->variation()->name eq $v_name, "pf 'variation' (using argument)");
ok($pf2->source()->name eq $source_name, "pf 'source' (using argument)");

# Test associated studies
my $pfs3 = $pf_adaptor->fetch_all_by_object_id('rs2299222');
my $asso_study = $pfs3->[0]->associated_studies->[0];
ok($asso_study->name eq 'asso_study', 'associated_studies');

# Test project_fullname
my $project = 'Full name of the project';
$pf2->project_fullname($project);
ok($pf2->project_fullname eq $project, 'project_fullname & _set_attribute');

# Test strain
ok($pf2->strain->name eq 'NA12891', 'strain');


done_testing();
