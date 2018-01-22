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
use Data::Dumper;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Test::MultiTestDB;
use Test::Exception;
our $verbose = 0;


my $omulti = Bio::EnsEMBL::Test::MultiTestDB->new('multi');
my $odb = $omulti->get_DBAdaptor('ontology');
Bio::EnsEMBL::Registry->add_db($omulti, 'ontology', $odb);

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');

my $vdb = $multi->get_DBAdaptor('variation');
my $db  = $multi->get_DBAdaptor('core');


$vdb->dnadb($db);

my $vfa = $vdb->get_VariationFeatureAdaptor();
my $va  = $vdb->get_VariationAdaptor();
my $vsa = $vdb->get_VariationSetAdaptor();
my $pa = $vdb->get_PopulationAdaptor();
my $source_adaptor = $vdb->get_SourceAdaptor(); 

ok($vfa && $vfa->isa('Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor'));

my $sa = $db->get_SliceAdaptor();

my $slice_phen    = $sa->fetch_by_region('chromosome','7');
my $slice_set     = $sa->fetch_by_region('chromosome','17');
my $slice_somatic = $sa->fetch_by_region('chromosome','13');
my $slice         = $sa->fetch_by_region('chromosome','18');

my $vfs = $vfa->fetch_all_by_Slice($slice);

my $vf_name  = 'rs142276873';
my $vf_id    = 33303674;


# Somatic
my $vf_somatic_name = 'COSM946275';

# Source
my $source_name = 'dbSNP';
my $somatic_source_name = 'COSMIC';
my $source = $source_adaptor->fetch_by_name($source_name);  
my $somatic_source = $source_adaptor->fetch_by_name($somatic_source_name); 

my $vf = $vfa->fetch_by_dbID($vf_id);

ok($vf->start() == 23821095,               "vf_id -> start");
ok($vf->end()   == 23821095,               "vf_id -> end") ;
ok($vf->strand() == 1,                     "vf_id -> strand");
ok($vf->allele_string() eq 'G/A',          "vf_id -> allele_string");
ok($vf->variation()->name() eq $vf_name,   "vf_id -> varname" );
ok($vf->display_id() eq $vf_name,          "vf_id -> display id");
ok($vf->map_weight() == 1,                 "vf_id -> map weight");
ok($vf->slice()->name() eq $slice->name(), "vf_id -> slice name");
ok($vf->display() ==1,                     "vf_id -> display=1");

print scalar @$vfs, "\n";
$vf = $vfs->[1];

ok($vf->dbID() == $vf_id,                  "var -> vf id");
ok($vf->slice->name() eq $slice->name(),   "var -> slice name ");
ok($vf->start == 23821095,                 "var -> start");
ok($vf->end() == 23821095,                 "var -> end");
ok($vf->strand() == 1,                     "var -> strand");
ok($vf->allele_string() eq 'G/A',          "var -> allele string");
ok($vf->variation()->name() eq $vf_name,   "var -> name");
ok($vf->display_id() eq $vf_name,          "var -> display id");
ok($vf->map_weight() == 1,                 "var -> map weight");
ok($vf->slice()->name() eq $slice->name(), "var -> slice name");
ok(@{$vf->consequence_type()} == 2 , "var -> consequence type count"); 
ok($vf->display_consequence() eq 'intron_variant', "var -> display consequence"); 
ok($vf->consequence_type()->[0] eq 'NMD_transcript_variant', "var -> consequence type"); 


my $cons = $vf->most_severe_OverlapConsequence();
ok($cons->description() eq 'A transcript variant occurring within an intron', "consequence description");
ok($cons->label() eq 'intron variant',                                        "consequence label");
ok($cons->SO_term() eq 'intron_variant',                                      "consequence SO_term"); 
ok($cons->SO_accession() eq 'SO:0001627',                                     "consequence SO_accession"); 
ok($cons->tier() eq '3',                                                      "consequence tier"); 
ok($cons->rank() eq '21',                                                     "consequence rank"); 
ok($cons->NCBI_term() eq 'intron',                                            "consequence NCBI term"); 
ok($cons->impact() eq 'MODIFIER',                                             "consequence impact"); 
ok($cons->display_term() eq 'INTRONIC',                                       "consequence display_term"); 


# test fetch_all_by_Variation inc failed
$va->db->include_failed_variations(1);
$vfa->db->include_failed_variations(1);

my $v = $va->fetch_by_dbID(14128071);
$vfs = $vfa->fetch_all_by_Variation($v);

ok(@$vfs == 1,                                 "var -> vf count ");
ok($vfs->[0]->display() eq '0',                "var -> display = 0");
ok($vfs->[0]->variation_name() eq 'rs67521280',"var -> vf_name" );
ok($vfs->[0]->dbID() eq 15275234 ,             "var -> vf_id" );


throws_ok { $vfa->fetch_all_by_Variation('Variant'); } qr/Variation arg expected/, 'fetch_all_by_Variation Throw on wrong argument';
throws_ok { $vfa->fetch_all_by_Variation(Bio::EnsEMBL::Variation::Variation->new()); } qr/Variation arg must have defined dbID/, 'fetch_all_by_Variation Throw on wrong argument';

my $vf2_name = 'rs2299222';

# test fetch all
print "\n# Test - fetch_all\n";
my $vfs2 = $vfa->fetch_all();
ok($vfs2->[0]->variation_name() eq $vf2_name, "vf by all");

# test fetch all somatic
print "\n# Test - fetch_all_somatic\n";
my $vfs3 = $vfa->fetch_all_somatic();
ok($vfs3->[0]->variation_name() eq $vf_somatic_name, "vf by all somatic");


## Slice ##

my $constraint = "vf.seq_region_start>20000000";

# test fetch all by Slice constraint with Variations
print "\n# Test - fetch_all_by_Slice_constraint_with_Variations\n";
my $vfs4 = $vfa->fetch_all_by_Slice_constraint_with_Variations($slice,$constraint);
ok($vfs4->[0]->variation_name() eq $vf_name, "vf by slice constraint with variation");

# test fetch all by Slice SO terms
print "\n# Test - fetch_all_by_Slice_SO_terms\n";
my $vfs4a = $vfa->fetch_all_by_Slice_SO_terms($slice);
ok(scalar @$vfs4a == 3, "fetch_all_by_Slice_SO_terms");
throws_ok { $vfa->fetch_all_by_Slice_SO_terms('Slice'); } qr/Slice arg expected/, 'fetch_all_by_Slice_SO_terms Throw on wrong argument';

$vfs4a = $vfa->fetch_all_by_Slice_SO_terms($slice, ['intron_variant']);
ok(scalar @$vfs4a == 2, "fetch_all_by_Slice_SO_terms");

# test fetch all by Slice constraint with TranscriptVariations
print "\n# Test - fetch_all_by_Slice_constraint_with_TranscriptVariations\n";
my $vfs5 = $vfa->fetch_all_by_Slice_constraint_with_TranscriptVariations($slice,$constraint);
ok($vfs5->[0]->variation_name() eq $vf_name, "vf by slice constraint with transcript variation");

# test fetch all genotyped by Slice
print "\n# Test - fetch_all_genotyped_by_Slice\n";
my $vfs6 = $vfa->fetch_all_genotyped_by_Slice($slice);
ok($vfs6->[0]->variation_name() eq $vf_name, "genotyped vf by slice");

## test fetch all with phenotype by Slice ##
print "\n# Test - fetch_all_with_phenotype_by_Slice\n";
my $vfs7 = $vfa->fetch_all_with_phenotype_by_Slice($slice_phen);
ok($vfs7->[0]->variation_name() eq $vf2_name, "vf with phenotype by slice");
throws_ok { $vfa->fetch_all_with_phenotype_by_Slice('Slice'); } qr/Slice arg expected/, 'fetch_all_with_phenotype_by_Slice Throw on wrong argument';

# fetch_all_with_phenotype_by_Slice - using variation source
my $vfs7a = $vfa->fetch_all_with_phenotype_by_Slice($slice_phen,'dbSNP');
ok($vfs7a->[0]->variation_name() eq $vf2_name, "vf with phenotype by slice - using variation source");

# fetch_all_with_phenotype_by_Slice - using phenotype source
my $vfs7b = $vfa->fetch_all_with_phenotype_by_Slice($slice_phen, undef, 'dbSNP');
ok($vfs7b->[0]->variation_name() eq $vf2_name, "vf with phenotype by slice - using phenotype source");

# fetch_all_with_phenotype_by_Slice - using phenotype
my $vfs7c = $vfa->fetch_all_with_phenotype_by_Slice($slice_phen, undef, undef, 'ACHONDROPLASIA');
ok($vfs7c->[0]->variation_name() eq $vf2_name, "vf with phenotype by slice - using phenotype name");
my $vfs7d = $vfa->fetch_all_with_phenotype_by_Slice($slice_phen, undef, undef, 1);
ok($vfs7d->[0]->variation_name() eq $vf2_name, "vf with phenotype by slice - using phenotype id");

## test fetch all with maf by Slice ##
print "\n# Test - fetch_all_with_maf_by_Slice\n";
# fetch_all_with_maf_by_Slice
my $vfs_maf_a = $vfa->fetch_all_with_maf_by_Slice($slice);
ok(scalar(@$vfs_maf_a) == 2, "vf with MAF by slice");

# fetch_all_with_maf_by_Slice - using a MAF threshold
my $maf = 0.1;
my $vfs_maf_b = $vfa->fetch_all_with_maf_by_Slice($slice,$maf);
ok($vfs_maf_b->[0]->minor_allele_frequency <= $maf, "vf with MAF lesser or equal than $maf by slice");
my $vfs_maf_c = $vfa->fetch_all_with_maf_by_Slice($slice,$maf,1);
ok($vfs_maf_c->[0]->minor_allele_frequency > $maf, "vf with MAF greater than $maf by slice");
throws_ok { $vfa->fetch_all_with_maf_by_Slice('Slice'); } qr/Slice arg expected/, 'fetch_all_with_maf_by_Slice Throw on wrong argument';
throws_ok { $vfa->fetch_all_with_maf_by_Slice($slice, 0.8); } qr/The maximum value of the minor allele frequency parameter should to be lesser/, 'fetch_all_with_maf_by_Slice MAF to big';

# test fetch all by Slice VariationSet
print "\n# Test - fetch_all_by_Slice_VariationSet\n";
my $vs = $vsa->fetch_by_name('1000 Genomes - All - common');
my $vfs8 = $vfa->fetch_all_by_Slice_VariationSet($slice_set,$vs);
ok($vfs8->[0]->variation_name() eq 'rs2255888', "vf by slice & variation set");
throws_ok { $vfa->fetch_all_by_Slice_VariationSet('Slice', $vs); } qr/Slice arg expected/, 'fetch_all_by_Slice_VariationSet Throw on wrong argument';
throws_ok { $vfa->fetch_all_by_Slice_VariationSet($slice_set, 'VariatioSet'); } qr/VariationSet arg expected/, 'fetch_all_by_Slice_VariationSet Throw on wrong argument';

# test fetch all by Slice Population
print "\n# Test - fetch_all_by_Slice_Population\n";
my $pop = $pa->fetch_by_name('SSMP:SSM');
my $vfs9 = $vfa->fetch_all_by_Slice_Population($slice,$pop);
ok($vfs9->[0]->variation_name() eq $vf_name, "vf by slice & population");
throws_ok { $vfa->fetch_all_by_Slice_Population('Slice', $pop); } qr/Slice arg expected/, 'fetch_all_by_Slice_Population Throw on wrong argument';
throws_ok { $vfa->fetch_all_by_Slice_Population($slice, 'Population'); } qr/Population arg expected/, 'fetch_all_by_Slice_Population Throw on wrong argument';

my $vfs9a = $vfa->fetch_all_by_Slice_Population($slice, $pop, 0.1);
ok($vfs9a->[0]->variation_name() eq $vf_name, "vf by slice & population & freq");

# test fetch all by Slice VariationSet & SO term
print "\n# Test - fetch_all_by_Slice_VariationSet & SO term\n";
my $vfs8b = $vfa->fetch_all_by_Slice_VariationSet_SO_terms($slice_set,$vs,['intron_variant'] );

ok($vfs8b->[0]->variation_name() eq 'rs182218163', "vf by slice & variation set & SO term");
throws_ok { $vfa->fetch_all_by_Slice_VariationSet_SO_terms('slice', $vs, ['intron_variant']); } qr/Slice arg expected/, 'Throw on wrong slice argument.';
throws_ok { $vfa->fetch_all_by_Slice_VariationSet_SO_terms($slice_set, 'variation_set', ['intron_variant']); } qr/VariationSet arg expected/, 'Throw on wrong variation set argument.';

my $vfs8c = $vfa->fetch_all_by_Slice_VariationSet_SO_terms($slice_set, $vs);
ok( scalar ( grep { $_->variation_name() eq 'rs2255888' } @$vfs8c) == 1, "vf by slice & variation set");

## Slice Somatic ##

# test fetch all somatic by Slice constraint
print "\n# Test - fetch_all_somatic_by_Slice_constraint\n";
my $vfs11 = $vfa->fetch_all_somatic_by_Slice_constraint($slice_somatic,$constraint);
ok($vfs11->[0]->variation_name() eq $vf_somatic_name, "somatic vf by slice constraint");

# test fetch all somatic by Slice
print "\n# Test - fetch_all_somatic_by_Slice\n";
my $vfs12 = $vfa->fetch_all_somatic_by_Slice($slice_somatic);
ok($vfs12->[0]->variation_name() eq $vf_somatic_name, "somatic vf by slice");

# fetch all somatic by Slice and Source
print "\n# Test - fetch_all_somatic_by_Slice_Source\n";
my $vfs12a = $vfa->fetch_all_somatic_by_Slice_Source($slice_somatic, $somatic_source);    
ok($vfs12a->[0]->variation_name() eq $vf_somatic_name, "somatic vf by slice and source");
throws_ok { $vfa->fetch_all_somatic_by_Slice_Source('slice', $somatic_source); } qr/Slice arg expected/, 'Throw on wrong slice argument.';
throws_ok { $vfa->fetch_all_somatic_by_Slice_Source($slice_somatic, 'somatic_source'); } qr/Source arg expected/, 'Throw on wrong source argument.';
my $somatic_source2 = Bio::EnsEMBL::Variation::Source->new(
  -name => 'test_source',
);
warns_like { $vfa->fetch_all_somatic_by_Slice_Source($slice_somatic, $somatic_source2); } qr/Source does not have dbID/, 'Warn on missing source dbID.';

# test fetch all somatic by Slice SO terms
print "\n# Test - fetch_all_somatic_by_Slice_SO_terms\n";
my $vfs12b = $vfa->fetch_all_somatic_by_Slice_SO_terms($slice_somatic, ['missense_variant']);
ok($vfs12b->[0]->variation_name() eq $vf_somatic_name, "somatic vf by slice and SO terms");
throws_ok { $vfa->fetch_all_somatic_by_Slice_SO_terms('slice', ['missense_variant']); } qr/Slice arg expected/, 'Throw on wrong slice argument.';
$vfs12b = $vfa->fetch_all_somatic_by_Slice_SO_terms($slice_somatic);
ok($vfs12b->[0]->variation_name() eq $vf_somatic_name, "somatic vf by slice and no SO terms");

# test fetch all somatic with phenotype by Slice
print "\n# Test - fetch_all_somatic_with_phenotype_by_Slice\n";
my $vfs13 = $vfa->fetch_all_somatic_with_phenotype_by_Slice($slice_somatic);
ok($vfs13->[0]->variation_name() eq $vf_somatic_name, "somatic vf with phenotype by slice");

# test fetch all somatic by Slice constraint with TranscriptVariations
print "\n# Test - fetch_all_somatic_by_Slice_constraint_with_TranscriptVariations\n";
my $vfs14 = $vfa->fetch_all_somatic_by_Slice_constraint_with_TranscriptVariations($slice_somatic,$constraint);
ok($vfs14->[0]->variation_name() eq $vf_somatic_name, "somatic vf by slice constraint with transcript variation");


## Other ##

# test fetch all with phenotype
print "\n# Test - fetch_all_with_phenotype\n";
my $vfs15 = $vfa->fetch_all_with_phenotype();
ok($vfs15->[0]->variation_name() eq $vf2_name, "vf with phenotype");

# fetch_all_with_phenotype - using variation source
my $vfs15a = $vfa->fetch_all_with_phenotype('dbSNP');
ok($vfs15a->[0]->variation_name() eq $vf2_name, "vf with phenotype - using variation source");

# fetch_all_with_phenotype - using phenotype source
my $vfs15b = $vfa->fetch_all_with_phenotype(undef, 'dbSNP');
ok($vfs15b->[0]->variation_name() eq $vf2_name, "vf with phenotype - using phenotype source");

# fetch_all_with_phenotype - using phenotype
my $vfs15c = $vfa->fetch_all_with_phenotype(undef, undef, 'ACHONDROPLASIA') ;
ok($vfs15c->[0]->variation_name() eq $vf2_name, "vf with phenotype - using phenotype");

# fetch_all_with_phenotype - using constraint
my $vfs15d = $vfa->fetch_all_with_phenotype(undef, undef, 'ACHONDROPLASIA', 'vf.seq_region_strand = 1') ;
ok($vfs15d->[0]->variation_name() eq $vf2_name, "vf with phenotype - using constraint");

# test fetch all somatic with phenotype
print "\n# Test - fetch_all_somatic_with_phenotype\n";
my $vfs16 = $vfa->fetch_all_somatic_with_phenotype();
ok($vfs16->[0]->variation_name() eq $vf_somatic_name, "somatic vf with phenotype");

# fetch_all_somatic_with_phenotype - using constraint
my $vfs16a = $vfa->fetch_all_somatic_with_phenotype(undef, undef, undef, 'vf.seq_region_strand = 1') ;
ok($vfs16a->[0]->variation_name() eq $vf_somatic_name, "somatic vf with phenotype - using constraint");

# test fetching VF with empty consequence type column
is($vfa->fetch_by_dbID(997738282)->display_consequence, 'sequence_variant', 'empty consequence column');


# test fetch Iterator
print "\n# Test - fetch_Iterator\n";
my $vfs17 = $vfa->fetch_Iterator();
ok($vfs17->next()->variation_name eq $vf2_name, 'vf fetch_Iterator');


# test list dbIDs
print "\n# Test - list_dbIDs\n";
my $dbIDs = $vfa->list_dbIDs();
ok(scalar(@$dbIDs) > 0 && $dbIDs->[0] == 1, 'vf list_dbIDs');

# test new fake
print "\n# Test - new_fake\n";
my $vfa_fake = Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor->new_fake('human');
ok($vfa_fake && $vfa_fake->isa('Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor'), 'fake vf adaptor - new_fake');

# store with missing _source_id
$source = Bio::EnsEMBL::Variation::Source->new(
  -name => 'dbSNP',
  -source_id => 1,
  -dbID => 1,
);
my $vf18 = Bio::EnsEMBL::Variation::VariationFeature->new(
  -start   => 100,
  -end     => 100,
  -strand  => 1,
  -slice   => $slice,
  -allele_string => 'A/T',
  -source => $source,
  -_variation_id => 1234,
  -is_somatic => 0,
  -class_attrib_id => 2,
  -class_SO_term => 'SNV',
);

$vfa->store($vf18);

# update with missing _source_id
my $upd_name = 'updated_test_vf18';
$vf18->variation_name($upd_name);
$vfa->update($vf18);
my $vf_new18 = $vfa->fetch_all_by_Slice_constraint($slice, "vf.variation_name='$upd_name'");
ok($vf_new18 && $vf_new18->[0]->variation_name eq $upd_name, "fetch updated vf");

# delete vf
my $dbID = $vf18->dbID;
my $dbh = $vfa->dbc->db_handle;
$dbh->do(qq{DELETE FROM variation_feature WHERE variation_feature_id=$dbID;}) or die $dbh->errstr;

## store ##
print "\n# Test - store\n";
delete $vf->{$_} for qw(dbID seq_region_start variation_name);
my $new_seq_region_start = 1000;
my $new_var_name = 'test';
$vf->start($new_seq_region_start);
$vf->variation_name("$new_var_name");

ok($vfa->store($vf), "store");

my $vfs_store = $vfa->fetch_all_by_Slice_constraint($slice, "vf.seq_region_start=$new_seq_region_start AND vf.variation_name='$new_var_name'");
ok($vfs_store && $vfs_store->[0]->seq_region_start == $new_seq_region_start && $vfs_store->[0]->variation_name eq $new_var_name, "fetch stored");

## update
print "\n# Test - update\n";
my $upd_vf = $vfs_store->[0];
$upd_name = 'updated_test';
$upd_vf->variation_name($upd_name);
$vfa->update($upd_vf);

my $vf_new = $vfa->fetch_all_by_Slice_constraint($slice, "vf.variation_name='$upd_name'");
ok($vf_new && $vf_new->[0]->variation_name eq $upd_name, "fetch updated vf");
$dbID = $vf_new->[0]->dbID;
$dbh->do(qq{DELETE FROM variation_feature WHERE variation_feature_id=$dbID;}) or die $dbh->errstr;

print "\n# Test - fetch_by_hgvs_notation\n";
my $hgvs_str = '9:g.139568335_1395683374GGCCGCTGGTGGGGATGGCTTCCAGCACCTGCACTGTGAC>GCGCAG';
throws_ok {$vfa->fetch_by_hgvs_notation($hgvs_str); } qr/Region requested must be smaller than 5kb/, 'Throw on region longer than 5kbt.';
done_testing();

