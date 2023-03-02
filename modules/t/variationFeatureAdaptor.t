# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2023] EMBL-European Bioinformatics Institute
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

# test fetch_all_by_Slice +/- inc failed with default caching
{
  print "\n# Test fetch all by Slice including failed flag, default cache\n";
  my $slice1 = $sa->fetch_by_region('chromosome','11',6303493,66324360);
  $vfa->db->include_failed_variations(1);
  my $vfs_slice = $vfa->fetch_all_by_Slice($slice1);
  cmp_ok(scalar @$vfs_slice, "==", 447, "slice (+failed default cache) -> vf count ");
  my $slice2 = $sa->fetch_by_region('chromosome','11',6303493,66324360);
  $vfa->db->include_failed_variations(0); # failed flag will be ignored as the cache is in place
  $vfs_slice = $vfa->fetch_all_by_Slice($slice2);
  cmp_ok(scalar @$vfs_slice,"==", 447, "slice (-failed default cache) -> vf count ");
}
# test fetch_all_by_Slice +/- inc failed with no caching
{
  print "\n# Test fetch all by Slice including failed flag, no caching\n";
  my $slice1 = $sa->fetch_by_region('chromosome','11',6303493,66324360);
  $vfa->db->include_failed_variations(1);
  $vfa->db->no_cache(1);
  my $vfs_slice = $vfa->fetch_all_by_Slice($slice1);
  cmp_ok(scalar @$vfs_slice, "==", 447, "slice (+failed, default cache) -> vf count ");
  my $slice2 = $sa->fetch_by_region('chromosome','11',6303493,66324360);
  $vfa->db->include_failed_variations(0); # failed flag will be included as the cache is not in place
  $vfs_slice = $vfa->fetch_all_by_Slice($slice2);
  cmp_ok(scalar @$vfs_slice,"==", 445, "slice (-failed, no cache) -> vf count ");
}
# test fetch_all_by_Slice +/- inc failed with clearing caching
{
  print "\n# Test fetch all by Slice including failed flag, clear cache\n";
  my $slice1 = $sa->fetch_by_region('chromosome','11',6303493,66324360);
  $vfa->db->include_failed_variations(1);
  $vfa->db()->no_cache(0);
  cmp_ok($vfa->db()->no_cache(), "==", 0, "caching is on");
  my $vfs_slice = $vfa->fetch_all_by_Slice($slice1); #will be saved in the cache
  cmp_ok(scalar @$vfs_slice, "==", 447, "slice (+failed, default cache) -> vf count ");
  my $slice2 = $sa->fetch_by_region('chromosome','11',6303493,66324360);
  $vfa->db->include_failed_variations(0); # failed flag will be ignored unless the cache is cleared
  $vfa->clear_cache(); # feature cache is cleared
  $vfs_slice = $vfa->fetch_all_by_Slice($slice2);
  cmp_ok(scalar @$vfs_slice,"==", 445, "slice (-failed, no cache) -> vf count ");
}

# test fetch all +/- inc failed
{
  print "\n# Test - fetch_all +/- inc failed\n";
  my $vf2_name = 'rs2299222';
  $vfa->db->include_failed_variations(0);
  my $vfs2 = $vfa->fetch_all();
  cmp_ok(scalar @$vfs2, "==", 1296, "vf by all - count (-failed)");
  cmp_ok($vfs2->[0]->variation_name(), "eq", $vf2_name, "vf by all - check first variation name");

  #test fetch all with inc failed my $vf_nameF='rs111067473';
  $vfa->db->include_failed_variations(1);
  my $vfs = $vfa->fetch_all(); 
  cmp_ok(scalar @$vfs, "==", 1303, "vf by all - count (+failed)");
}

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

my $vfs8d = $vfa->fetch_all_by_Slice_Source($slice_set, $source);
ok(scalar @$vfs8d == 3, "fetch_all_by_Slice_Source");
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

ok(scalar @{$vfa->fetch_all_by_location_identifier('18:40228819:A_G')} == 1, "fetch_all_by_location_identifier '18:40228819:A_G'");

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
  -adaptor => $vfa,
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
$hgvs_str = 'Q00872:p.Ala53Val';
ok($vfa->fetch_by_hgvs_notation($hgvs_str)->allele_string eq 'C/T', 'HGVSp notation using UniProt ID');
ok($vfa->fetch_by_hgvs_notation('ENST00000470094:c.55_111del')->end eq 32954180, 'HGVSc multi-exon deletion');

# test HGVS protein when codon is within two exons
# forward strand
ok($vfa->fetch_by_hgvs_notation('ENSP00000422007.1:p.Gly469Glu')->start eq 66325646, 'HGVSp multi-exon (forward)');
# reverse strand
ok($vfa->fetch_by_hgvs_notation('ENSP00000293261:p.Arg232Met')->start eq 48846578, 'HGVSp multi-exon (reverse)');


print "\n# Test - fetch_by_spdi_notation\n";
my $spdi_str = 'NC_000013.10:32954017::';
throws_ok {$vfa->fetch_by_spdi_notation($spdi_str); } qr/Could not parse the SPDI notation $spdi_str/, 'Throw on invalid SPDI notation.';
$spdi_str = 'N_000013.10:32954017:G:A';
throws_ok {$vfa->fetch_by_spdi_notation($spdi_str); } qr/Sequence name N_000013.10 not valid/, 'Throw on invalid sequence id.';
$spdi_str = 'NC_000013.10:32954017:G:A:';
throws_ok {$vfa->fetch_by_spdi_notation($spdi_str); } qr/Could not parse the SPDI notation $spdi_str. Too many elements present/, 'Throw on invalid SPDI notation. Too many elements.'; 
$spdi_str = 'NC_000013.10:32954017:G';
throws_ok {$vfa->fetch_by_spdi_notation($spdi_str); } qr/Could not parse the SPDI notation $spdi_str. Too few elements present/, 'Throw on invalid SPDI notation. Too few elements.'; 
$spdi_str = 'NC_000013:32954017:C:A';
throws_ok {$vfa->fetch_by_spdi_notation($spdi_str); } qr/Reference allele extracted from NC_000013:32954018-32954018 \(G\) does not match reference allele given by SPDI notation $spdi_str \(C\)/, 'Throw on reference allele does not match SPDI reference allele.'; 
$spdi_str = 'NC_000013:32954017:2:A';
throws_ok {$vfa->fetch_by_spdi_notation($spdi_str); } qr/Could not parse the SPDI notation $spdi_str. Deleted sequence length \(2\) does not match inserted sequence length \(1\)./, 'Throw on invalid SPDI notation. Deleted sequence wrong size.'; 
$spdi_str = 'NC_000013.10:32954017:G:G';
throws_ok {$vfa->fetch_by_spdi_notation($spdi_str); } qr/Reference allele given by SPDI notation $spdi_str \(G\) matches alt allele given by SPDI notation $spdi_str \(G\)/, 'Throw on invalid alt allele.'; 
$spdi_str = 'NC_000011.9:66317226::1';
throws_ok {$vfa->fetch_by_spdi_notation($spdi_str); } qr/Could not parse the SPDI notation $spdi_str. SPDI notation not supported./, 'Throw on invalid SPDI notation. Insertion format not valid.'; 
$spdi_str = 'NC_000013:32954017:0:0';
throws_ok {$vfa->fetch_by_spdi_notation($spdi_str); } qr/Could not parse the SPDI notation $spdi_str. SPDI notation not supported./, 'Throw on invalid SPDI notation.'; 
$spdi_str = 'NC_000013:32954017:G1:0';
throws_ok {$vfa->fetch_by_spdi_notation($spdi_str); } qr/Could not parse the SPDI notation $spdi_str. SPDI notation not supported./, 'Throw on invalid SPDI notation. Deleted sequence has numbers and digits.'; 
$spdi_str = 'NC_000013.10:32954017:G:A';
$vf = $vfa->fetch_by_spdi_notation($spdi_str);
ok($vf->seq_region_start eq '32954018', "Valid substitution 'NC_000013.10:32954017:G:A', spdi position is 0-based: 32954017(spdi) = 32954018");
$spdi_str = 'NC_000013.10:32954017:1:A';
$vf = $vfa->fetch_by_spdi_notation($spdi_str);
ok($vf->seq_region_start eq '32954018', "Valid substitution 'NC_000013.10:32954017:1:A', spdi position is 0-based: 32954017(spdi) = 32954018");
$spdi_str = 'NC_000011.9:66321302:TG:CA';
$vf = $vfa->fetch_by_spdi_notation($spdi_str);
ok($vf->seq_region_start eq '66321303' && $vf->seq_region_end eq '66321304', "Valid substitution 'NC_000011.9:66321302:TG:CA'");
$spdi_str = 'NC_000002.11:45406939:N:';
$vf = $vfa->fetch_by_spdi_notation($spdi_str);
ok($vf->allele_string eq 'N/-', "Valid deletion 'NC_000002.11:45406939:N:' - reference and spdi alleles match"); 
$spdi_str = 'NC_000002.11:45406939:1:';
$vf = $vfa->fetch_by_spdi_notation($spdi_str);
ok($vf->allele_string eq 'N/-', "Valid deletion 'NC_000002.11:45406939:1:' - reference and spdi alleles match");
$spdi_str = 'NC_000013:32954017:G:0';
$vf = $vfa->fetch_by_spdi_notation($spdi_str); 
ok($vf->allele_string eq 'G/-', "Valid deletion 'NC_000013:32954017:G:0'");
$spdi_str = 'NC_000011:66320351:10:0';
$vf = $vfa->fetch_by_spdi_notation($spdi_str); 
ok($vf->allele_string eq 'CACACACACA/-', "Valid deletion 'NC_000011:66320351:10:0'");
$spdi_str = 'NC_000012:101997654:GCCATGATCATGCCACTGCACTCCT:';
$vf = $vfa->fetch_by_spdi_notation($spdi_str); 
ok($vf->allele_string eq 'GCCATGATCATGCCACTGCACTCCT/-', "Valid deletion 'NC_000012:101603876:GCCATGATCATGCCACTGCACTCCT:'");  
$spdi_str = 'NC_000011.9:66317226::C';
$vf = $vfa->fetch_by_spdi_notation($spdi_str);
ok($vf->seq_region_start eq '66317227' && $vf->seq_region_end eq '66317226', "Valid insertion 'NC_000011.9:66317226::C' - a 'C' is inserted at 66317226-66317227");
$spdi_str = 'NC_000013:32954017:0:A';
$vf = $vfa->fetch_by_spdi_notation($spdi_str); 
ok($vf->allele_string eq '-/A', "Valid insertion 'NC_000013:32954017:0:A'"); 
$spdi_str = 'NC_000011.9:66321302:2:CA';
$vf = $vfa->fetch_by_spdi_notation($spdi_str);
ok($vf->seq_region_start eq '66321303' && $vf->seq_region_end eq '66321304', "Valid substitution 'NC_000011.9:66321302:2:CA'"); 
$spdi_str = 'NC_000012:102009450:GCCCCCC:CT'; 
$vf = $vfa->fetch_by_spdi_notation($spdi_str);
ok($vf->allele_string eq 'GCCCCCC/CT' , "Valid indel"); 
my $spdi_lrg = 'LRG_293:69534:N:G';
$vf = $vfa->fetch_by_spdi_notation($spdi_lrg);
ok($vf->allele_string eq 'N/G', "LRG - Valid insertion 'LRG_293:69534:N:G'");
$vf = $vfa->fetch_by_spdi_notation('NT_004487.20:127830:N:C');
ok($vf->allele_string eq 'N/C', "NT - Valid substitution 'NT_004487.20:127830:N:C'");

## check ref matching
my $bad_hgvs = 'ENSP00000434898.1:p.Cys6Ser';
throws_ok {$vfa->fetch_by_hgvs_notation($bad_hgvs); }qr/Sequence translated from reference \(TCT -> S\) does not match input sequence \(C\)/, 'Throw if HGVS does not match reference';

throws_ok {$vfa->fetch_by_hgvs_notation('ENST00000470094.1:c.59A>G'); }qr/Reference allele extracted from ENST00000470094:32954035-32954035 \(G\) does not match reference allele given by HGVS notation ENST00000470094\.1:c\.59A>G \(A\)/, 'Throws if reference allele does not match';

my $ok_hgvs = 'ENSP00000293261.2:p.Ser455del';
$vf = $vfa->fetch_by_hgvs_notation($ok_hgvs);
ok($vf->allele_string eq 'AGC/-', "HGVSp matches reference");

my $lrg_hgvsg = 'LRG_293:g.69535N>G';
$vf = $vfa->fetch_by_hgvs_notation($lrg_hgvsg);
ok($vf->allele_string eq 'N/G', "HGVSg LRG");

print "\n# Test - fetch_all_by_caid\n";
my $caid;

throws_ok {$vfa->fetch_all_by_caid(); } qr/No CAid provided/,
  'Throw on no CAid';

throws_ok {$vfa->fetch_all_by_caid(''); } qr/No CAid provided/,
  'Throw on empty CAid';

throws_ok {$vfa->fetch_all_by_caid('CAXXX'); } qr/CAid has an invalid format/,
  'Throw on invalid CAid format';

$caid = 'CA124';
$vfs = $vfa->fetch_all_by_caid($caid);
ok(@$vfs == 0, "CAid ($caid) - no vf found");

# CAid for biallelic SNV
$caid = 'CA6124251';
$vfs = $vfa->fetch_all_by_caid($caid);
ok(@$vfs == 1, "CAid ($caid) - biallelic single vf found");
cmp_ok($vfs->[0]->variation_name(), 'eq', 'rs490998', "CAid ($caid) - biallelic single vf -> vf_name");
cmp_ok($vfs->[0]->allele_string(), 'eq', 'G/A', "CAid ($caid) - biallelic single vf -> vf_allele_string");

# CAid for multi-allelic SNV
$caid = 'CA6124851';
$vfs = $vfa->fetch_all_by_caid($caid);
ok(@$vfs == 1, "CAid ($caid) - multi-allelic single vf found");
cmp_ok($vfs->[0]->variation_name(), 'eq', 'rs377076795', "CAid ($caid) - multi-allelic single vf -> vf_name");
cmp_ok($vfs->[0]->allele_string(), 'eq', 'C/T', "CAid ($caid) - multi-allelic single vf -> vf_allele_string");

# CAid for multiple variants
$caid = 'CA169452307';
$vfs = $vfa->fetch_all_by_caid($caid);
my @vfs_sort = sort {$a->variation_name() cmp $b->variation_name} @{$vfs};
ok(@vfs_sort == 2, "CAid ($caid) - multiple vf found");
cmp_ok($vfs_sort[0]->variation_name(), 'eq', 'rs7811371', "CAid ($caid) - first vf -> vf_name");
cmp_ok($vfs_sort[1]->variation_name(), 'eq', 'rs7811371_dup', "CAid ($caid) - second vf -> vf_name");
cmp_ok($vfs_sort[0]->allele_string(), 'eq', 'G/A', "CAid ($caid) - first vf -> vf_allele_string");
cmp_ok($vfs_sort[1]->allele_string(), 'eq', 'G/A', "CAid ($caid) - second vf -> vf_allele_string");

# CAid for failed variant
$caid = 'CA025949';
$vfa->db->include_failed_variations(0);
$caid = 'CA025949';
$vfs = $vfa->fetch_all_by_caid($caid);
ok(@$vfs == 0, "CAid ($caid) - exclude failed vf");

$vfa->db->include_failed_variations(1);
$caid = 'CA025949';
$vfs = $vfa->fetch_all_by_caid($caid);
ok(@$vfs == 1, "CAid ($caid) - include failed vf");
cmp_ok($vfs->[0]->variation_name(), 'eq', 'rs80359157', "CAid ($caid) - failed vf -> vf_name");
cmp_ok($vfs->[0]->allele_string(), 'eq', 'C/G', "CAid ($caid) - failed vf -> vf_allele_string");

done_testing();
