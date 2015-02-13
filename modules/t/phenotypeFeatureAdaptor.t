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

# fetch_all_by_object_id
my $pfs = $pfa->fetch_all_by_object_id('rs2299222');
ok(ref($pfs) eq 'ARRAY' && scalar @$pfs == 1 && $pfs->[0]->object_id eq 'rs2299222', "fetch_all_by_object_id");

# fetch_all_by_Slice_type
my $sla = $multi->get_DBAdaptor('core')->get_SliceAdaptor;
my $sl  = $sla->fetch_by_region('chromosome', 7, 86442403, 86442405);
$pfs = $pfa->fetch_all_by_Slice_type($sl, 'Variation');
ok(ref($pfs) eq 'ARRAY' && scalar @$pfs == 1 && $pfs->[0]->object_id eq 'rs2299222', "fetch_all_by_Slice_type");

# fetch_all_by_Variation
my $va = $vdba->get_VariationAdaptor();
my $v  = $va->fetch_by_name('rs2299222');
$pfs = $pfa->fetch_all_by_Variation($v);
ok(ref($pfs) eq 'ARRAY' && scalar @$pfs == 1 && $pfs->[0]->object_id eq 'rs2299222', "fetch_all_by_Variation");

# fetch_all_by_Variation_list
$pfs = $pfa->fetch_all_by_Variation_list([$v]);
ok(ref($pfs) eq 'ARRAY' && scalar @$pfs == 1 && $pfs->[0]->object_id eq 'rs2299222', "fetch_all_by_Variation_list");

# fetch_all_by_StructuralVariation
my $sva = $vdba->get_StructuralVariationAdaptor();
my $sv  = $sva->fetch_by_name('esv2751608');
$pfs = $pfa->fetch_all_by_StructuralVariation($sv);
ok(ref($pfs) eq 'ARRAY' && scalar @$pfs == 1 && $pfs->[0]->object_id eq 'esv2751608', "fetch_all_by_StructuralVariation");

# fetch_all_by_Gene
my $ga = $multi->get_DBAdaptor('core')->get_GeneAdaptor;
my $g  = $ga->fetch_by_stable_id('ENSG00000176105');
$pfs = $pfa->fetch_all_by_Gene($g);
ok(ref($pfs) eq 'ARRAY' && scalar @$pfs == 1 && $pfs->[0]->object_id eq 'ENSG00000176105', "fetch_all_by_Gene");

# fetch_all_by_VariationFeature_list
$pfs = $pfa->fetch_all_by_VariationFeature_list($v->get_all_VariationFeatures);
ok(ref($pfs) eq 'ARRAY' && scalar @$pfs == 1 && $pfs->[0]->object_id eq 'rs2299222', "fetch_all_by_VariationFeature_list");

# fetch_all_by_Study
my $sa = $vdba->get_StudyAdaptor();
my $s  = $sa->fetch_by_dbID(4237);
$pfs = $pfa->fetch_all_by_Study($s);
ok(ref($pfs) eq 'ARRAY' && scalar @$pfs == 2 && (grep {$_->object_id eq 'rs2299222'} @$pfs), "fetch_all_by_Study");

# fetch_all_by_phenotype_name_source_name
$pfs = $pfa->fetch_all_by_phenotype_name_source_name('ACH', 'dbSNP');
ok(ref($pfs) eq 'ARRAY' && scalar @$pfs == 2 && (grep {$_->object_id eq 'rs2299222'} @$pfs), "fetch_all_by_phenotype_name_source_name");

# fetch_all_by_phenotype_description_source_name
$pfs = $pfa->fetch_all_by_phenotype_description_source_name('ACHONDROPLASIA', 'dbSNP');
ok(ref($pfs) eq 'ARRAY' && scalar @$pfs == 2 && (grep {$_->object_id eq 'rs2299222'} @$pfs), "fetch_all_by_phenotype_description_source_name");

# fetch_all_by_phenotype_id_source_name
$pfs = $pfa->fetch_all_by_phenotype_id_source_name(1, 'dbSNP');
ok(ref($pfs) eq 'ARRAY' && scalar @$pfs == 2 && (grep {$_->object_id eq 'rs2299222'} @$pfs), "fetch_all_by_phenotype_id_source_name");

# fetch_all_by_associated_gene
$pfs = $pfa->fetch_all_by_associated_gene('YES1');
ok(ref($pfs) eq 'ARRAY' && scalar @$pfs == 1 && (grep {$_->object_id eq 'rs2299222'} @$pfs), "fetch_all_by_associated_gene");

# count_all_by_associated_gene
my $count = $pfa->count_all_by_associated_gene('YES1');
ok($count && $count == 1, "count_all_by_associated_gene");

# count_all_by_Phenotype
my $pa = $vdba->get_PhenotypeAdaptor();
my $p  = $pa->fetch_by_dbID(1);
$count = $pfa->count_all_by_Phenotype($p);
ok($count && $count == 4, "count_all_by_Phenotype");

# count_all_by_Gene
$count = $pfa->count_all_by_Gene($g);
ok($count && $count == 1, "count_all_by_Gene");

# count_all_by_phenotype_id
$count = $pfa->count_all_by_phenotype_id(1);
ok($count && $count == 4, "count_all_by_phenotype_id");

# fetch_all
$pfs = $pfa->fetch_all();
ok(ref($pfs) eq 'ARRAY' && scalar @$pfs == 5 && (grep {$_->object_id eq 'rs2299222'} @$pfs), "fetch_all");

# store
my $pf = $pfs->[0];
delete $pf->{dbID};
$pf->object_id('test');
$pf->{source_id} = 1;
ok($pfa->store($pf), "store");

$pfs = $pfa->fetch_all_by_object_id('test');
ok(ref($pfs) eq 'ARRAY' && scalar @$pfs == 1 && (grep {$_->object_id eq 'test'} @$pfs), "fetch stored");

done_testing();

