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

use Test::Exception;
use Test::More;
use FindBin qw($Bin);

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::Population;
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Test::MultiTestDB;
use FileHandle;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');
my $vdba = $multi->get_DBAdaptor('variation');

my $pfa = $vdba->get_PhenotypeFeatureAdaptor();
ok($pfa && $pfa->isa('Bio::EnsEMBL::Variation::DBSQL::PhenotypeFeatureAdaptor'), "get PhenotypeFeature adaptor");

# Slice
my $sla = $multi->get_DBAdaptor('core')->get_SliceAdaptor;
my $sl  = $sla->fetch_by_region('chromosome', 7, 86442403, 86442405);

# fetch_all_by_object_id
my $pfs = $pfa->fetch_all_by_object_id('rs2299222');
ok(ref($pfs) eq 'ARRAY' && scalar @$pfs == 1 && $pfs->[0]->object_id eq 'rs2299222', "fetch_all_by_object_id");
$pfs =  $pfa->fetch_all_by_object_id('rs2299222', 'Variation');
ok(ref($pfs) eq 'ARRAY' && scalar @$pfs == 1 && $pfs->[0]->object_id eq 'rs2299222', "fetch_all_by_object_id + type");
throws_ok { $pfa->fetch_all_by_object_id('rs2299222', 'Variant'); } qr/is not a valid object type, valid types are/, ' > Throw on wrong object type';

# fetch_all_by_Slice_type
$pfs = $pfa->fetch_all_by_Slice_type($sl, 'Variation');
ok(ref($pfs) eq 'ARRAY' && scalar @$pfs == 1 && $pfs->[0]->object_id eq 'rs2299222', "fetch_all_by_Slice_type");

# fetch_all_by_Variation
my $va = $vdba->get_VariationAdaptor();
my $v  = $va->fetch_by_name('rs2299222');
$pfs = $pfa->fetch_all_by_Variation($v);
ok(ref($pfs) eq 'ARRAY' && scalar @$pfs == 1 && $pfs->[0]->object_id eq 'rs2299222', "fetch_all_by_Variation");
throws_ok { $pfa->fetch_all_by_Variation('Variant'); } qr/Variation arg expected/, ' > Throw on wrong argument';


# fetch_all_by_Variation_list
$pfs = $pfa->fetch_all_by_Variation_list([$v]);
ok(ref($pfs) eq 'ARRAY' && scalar @$pfs == 1 && $pfs->[0]->object_id eq 'rs2299222', "fetch_all_by_Variation_list");
throws_ok { $pfa->fetch_all_by_Variation_list(['Variant']); } qr/Variation arg expected/, ' > Throw on wrong argument';
throws_ok { $pfa->fetch_all_by_Variation_list([Bio::EnsEMBL::Variation::Variation->new()]); } qr/Variation arg must have defined name/, ' > Throw on wrong argument';

# fetch_all_by_StructuralVariation
my $sva = $vdba->get_StructuralVariationAdaptor();
my $sv  = $sva->fetch_by_name('esv2751608');
$pfs = $pfa->fetch_all_by_StructuralVariation($sv);
ok(ref($pfs) eq 'ARRAY' && scalar @$pfs == 1 && $pfs->[0]->object_id eq 'esv2751608', "fetch_all_by_StructuralVariation");
throws_ok { $pfa->fetch_all_by_StructuralVariation('Variant'); } qr/BaseStructuralVariation arg expected/, ' > Throw on wrong argument';

# fetch_all_by_Gene
my $ga = $multi->get_DBAdaptor('core')->get_GeneAdaptor;
my $g  = $ga->fetch_by_stable_id('ENSG00000176105');
$pfs = $pfa->fetch_all_by_Gene($g);
ok(ref($pfs) eq 'ARRAY' && scalar @$pfs == 1 && $pfs->[0]->object_id eq 'ENSG00000176105', "fetch_all_by_Gene");
throws_ok { $pfa->fetch_all_by_Gene('gene'); } qr/Gene arg expected/, ' > Throw on wrong argument';

# fetch_all_by_VariationFeature_list
$pfs = $pfa->fetch_all_by_VariationFeature_list($v->get_all_VariationFeatures);
ok(ref($pfs) eq 'ARRAY' && scalar @$pfs == 1 && $pfs->[0]->object_id eq 'rs2299222', "fetch_all_by_VariationFeature_list");
throws_ok { $pfa->fetch_all_by_VariationFeature_list(['VariationFeature']); } qr/VariationFeature arg expected/, ' > Throw on wrong argument';
throws_ok { $pfa->fetch_all_by_VariationFeature_list([Bio::EnsEMBL::Variation::VariationFeature->new()]); } qr/VariationFeatures in list must have defined names/, ' > Throw on wrong argument';

# fetch_all_by_Study
my $sa = $vdba->get_StudyAdaptor();
my $s  = $sa->fetch_by_dbID(4237);
$pfs = $pfa->fetch_all_by_Study($s);
ok(ref($pfs) eq 'ARRAY' && scalar @$pfs == 2 && (grep {$_->object_id eq 'rs2299222'} @$pfs), "fetch_all_by_Study");
throws_ok { $pfa->fetch_all_by_Study('Study'); } qr/Study arg expected/, ' > Throw on wrong argument';
throws_ok { $pfa->fetch_all_by_Study(Bio::EnsEMBL::Variation::Study->new()); } qr/Study arg must have defined dbID/, ' > Throw on wrong argument';

# fetch_all_by_Slice_Study
$pfs = $pfa->fetch_all_by_Slice_Study($sl,$s);
ok(ref($pfs) eq 'ARRAY' && scalar @$pfs == 1 && (grep {$_->object_id eq 'rs2299222'} @$pfs), "fetch_all_by_Slice_Study");
throws_ok { $pfa->fetch_all_by_Slice_Study('Slice', $s); } qr/Slice arg expected/, ' > Throw on wrong argument';
throws_ok { $pfa->fetch_all_by_Slice_Study($sl,'Study'); } qr/Study arg expected/, ' > Throw on wrong argument';
throws_ok { $pfa->fetch_all_by_Slice_Study($sl,Bio::EnsEMBL::Variation::Study->new()); } qr/Study arg must have defined dbID/, ' > Throw on wrong argument';

# fetch_all_by_phenotype_name_source_name
$pfs = $pfa->fetch_all_by_phenotype_name_source_name('ACH', 'dbSNP');
ok(ref($pfs) eq 'ARRAY' && scalar @$pfs == 2 && (grep {$_->object_id eq 'rs2299222'} @$pfs), "fetch_all_by_phenotype_name_source_name");

# fetch_all_by_phenotype_description_source_name
$pfs = $pfa->fetch_all_by_phenotype_description_source_name('ACHONDROPLASIA', 'dbSNP');
ok(ref($pfs) eq 'ARRAY' && scalar @$pfs == 2 && (grep {$_->object_id eq 'rs2299222'} @$pfs), "fetch_all_by_phenotype_description_source_name");

# fetch_all_by_phenotype_id_source_name
$pfs = $pfa->fetch_all_by_phenotype_id_source_name(1, 'dbSNP');
ok(ref($pfs) eq 'ARRAY' && scalar @$pfs == 2 && (grep {$_->object_id eq 'rs2299222'} @$pfs), "fetch_all_by_phenotype_id_source_name");

# fetch_all_by_phenotype_id_feature_type
$pfs = $pfa->fetch_all_by_phenotype_id_feature_type(1, 'Gene');
ok(ref($pfs) eq 'ARRAY' && scalar @$pfs == 1 && (grep {$_->object_id eq 'ENSG00000176105'} @$pfs), "fetch_all_by_phenotype_id_feature_type");

# fetch_all_by_Slice_with_ontology_accession
my $sl_oa  = $sla->fetch_by_region('chromosome', 13, 86442400, 86442450);
$pfs = $pfa->fetch_all_by_Slice_with_ontology_accession($sl_oa, 'Variation');

ok(ref($pfs) eq 'ARRAY' && scalar @$pfs == 1 &&  $pfs->[0]->object_id eq 'rs2299299' && $pfs->[0]->get_all_ontology_accessions->[0] eq 'Orphanet:130', "fetch_all_by_Slice_with_ontology_accession");

# fetch_all_by_phenotype_ontology_accession
$pfs = $pfa->fetch_all_by_phenotype_accession_source('Orphanet:130');
ok(ref($pfs) eq 'ARRAY' && scalar @$pfs == 1 && (grep {$_->object_id eq 'rs2299299'} @$pfs), "fetch_all_by_phenotype_accession");

# fetch_all_by_phenotype_ontology_accession + source
$pfs = $pfa->fetch_all_by_phenotype_accession_source('Orphanet:130', 'dbSNP');
ok(ref($pfs) eq 'ARRAY' && scalar @$pfs == 0 , "fetch_all_by_phenotype_accession + source");

# fetch_all_by_phenotype_ontology_accession & map type
$pfs = $pfa->fetch_all_by_phenotype_accession_type_source('Orphanet:130', 'is');
ok(ref($pfs) eq 'ARRAY' && scalar @$pfs == 1 && (grep {$_->object_id eq 'rs2299299'} @$pfs), "fetch_all_by_phenotype_accession_type_source");

# fetch_all_by_phenotype_ontology_accession & map type
$pfs = $pfa->fetch_all_by_phenotype_accession_type_source('Orphanet:130', 'involves');
ok(ref($pfs) eq 'ARRAY' && scalar @$pfs == 0 , "fetch_all_by_phenotype_accession_type_source");


# fetch_all_by_associated_gene
$pfs = $pfa->fetch_all_by_associated_gene('YES1');
ok(ref($pfs) eq 'ARRAY' && scalar @$pfs == 1 && (grep {$_->object_id eq 'rs2299222'} @$pfs), "fetch_all_by_associated_gene");

# fetch_all_by_associated_gene_phenotype_description
$pfs = $pfa->fetch_all_by_associated_gene_phenotype_description('YES1', 'ACHONDROPLASIA');
ok(ref($pfs) eq 'ARRAY' && scalar @$pfs == 1 && (grep {$_->object_id eq 'rs2299222'} @$pfs), "fetch_all_by_associated_gene_phenotype_description");

# fetch_all_by_Phenotype
my $pa = $vdba->get_PhenotypeAdaptor();
my $p  = $pa->fetch_by_dbID(1);
$pfs = $pfa->fetch_all_by_Phenotype($p);
ok($pfs->[0]->object_id() eq 'rs2299222', "fetch_all_by_Phenotype") ;
throws_ok { $pfa->fetch_all_by_Phenotype(); } qr/Phenotype arg expected/, ' > Throw on missing argument';

# count_all_by_associated_gene
my $count = $pfa->count_all_by_associated_gene('YES1');
ok($count && $count == 1, "count_all_by_associated_gene");

# count_all_by_Phenotype
$count = $pfa->count_all_by_Phenotype($p);
ok($count && $count == 4, "count_all_by_Phenotype");

# count_all_by_Gene
$count = $pfa->count_all_by_Gene($g);
ok($count && $count == 1, "count_all_by_Gene");

# count_all_by_phenotype_id
$count = $pfa->count_all_by_phenotype_id(1);
ok($count && $count == 4, "count_all_by_phenotype_id");

# count_all_type_by_phenotype_id
$count = $pfa->count_all_type_by_phenotype_id(1);
ok($count && $count->{'Variation'} == 1 && $count->{'StructuralVariation'} == 2 && $count->{'Gene'} == 1, "count_all_type_by_phenotype_id");

# fetch_all
$pfs = $pfa->fetch_all();
ok(ref($pfs) eq 'ARRAY' && scalar @$pfs == 6 && (grep {$_->object_id eq 'rs2299222'} @$pfs), "fetch_all");

# store
my $pf = $pfs->[0];
delete $pf->{dbID};
$pf->object_id('test');
$pf->{source_id} = 1;
ok($pfa->store($pf), "store");

$pfs = $pfa->fetch_all_by_object_id('test');
ok(ref($pfs) eq 'ARRAY' && scalar @$pfs == 1 && (grep {$_->object_id eq 'test'} @$pfs), "fetch stored");

$count = $pfa->_check_gene_by_HGNC('YES1'); 
ok($count == 1, '_check_gene_by_HGNC');

my $attribs = $pfa->_fetch_attribs_by_dbID(1);
ok($attribs->{'associated_gene'} eq 'YES1', '_fetch_attribs_by_dbID');
throws_ok { $pfa->_fetch_attribs_by_dbID } qr/Cannot fetch attributes without dbID/, ' > Throw on missing dbID';


### test trailing white space removal on attrib value
my $padded_genename  = "gene name  ";
my $clipped_genename = "gene name";
my $v_name = 'rs_ws';
my $submitter_name = 'lab name';
my $last_eval_date = '2010-10-10';

my $variation = Bio::EnsEMBL::Variation::Variation->new(-name   => $v_name,
                                                        -source => $pf->source);

my $pf_ws = Bio::EnsEMBL::Variation::PhenotypeFeature->new(
    -slice     => $sl,
    -start     => 23821095,
    -end       => 23821095,
    -phenotype => $pf->phenotype(),
    -type      => 'Variation',
    -object    => $variation,
    -source    => $pf->source(),
    -is_significant =>1,
    -attribs   => {
      associated_gene   => $padded_genename,
      submitter_names   => [$submitter_name],
      DateLastEvaluated => $last_eval_date
    },
);


$pfa->store($pf_ws);
my $pf_ws_ext = $pfa->fetch_all_by_object_id('rs_ws');
ok($pf_ws_ext->[0]->associated_gene() eq $clipped_genename, "trailing white space removal" );
ok($pf_ws_ext->[0]->submitter_names()->[0] eq $submitter_name, "submitter name stored & retrieved");
ok($pf_ws_ext->[0]->date_last_evaluated() eq $last_eval_date, "last evaluation date stored & retrieved");

done_testing();

