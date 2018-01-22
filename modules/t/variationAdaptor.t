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

use Bio::EnsEMBL::Variation::DBSQL::VariationSetAdaptor;
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Test::MultiTestDB;


our $verbose = 0;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');

my $vdb = $multi->get_DBAdaptor('variation');

my $va = $vdb->get_VariationAdaptor();

ok($va && $va->isa('Bio::EnsEMBL::Variation::DBSQL::VariationAdaptor'), "isa var adaptor");

my $vsa = $vdb->get_VariationSetAdaptor();
my $pa  = $vdb->get_PopulationAdaptor();

my $var_name = 'rs7569578';
my $var_id   = 4770800;

# test fetch by dbID

my $var = $va->fetch_by_dbID($var_id);

ok($var->name() eq $var_name,      'variation name by var id');
ok($var->stable_id() eq $var_name,      'stable id by var id');
ok($var->source_name() eq 'dbSNP',          'source name by var id'   );
ok($var->source_description() eq 'Variants (including SNPs and indels) imported from dbSNP',   'source description by var id'   );
ok($var->source_version() == 138, 'source version by var id');
ok($var->source_url() eq 'http://www.ncbi.nlm.nih.gov/projects/SNP/', 'source url by var id');
ok($va->get_source_version('dbSNP') == 138, 'source version by var adaptor/name');
ok($var->ancestral_allele() eq 'A',  'ancestral_allele by var id');
ok($var->minor_allele() eq 'A',      'minor allele by var id');
ok($var->minor_allele_count() == 358, 'minor allele count by var id');
ok($var->minor_allele_frequency() eq '0.164371', 'minor allele frequency by var id' );
ok($var->get_all_clinical_significance_states()->[0] eq 'benign', 'clinsig by var id');
ok($var->display_consequence() eq 'intron_variant', 'display SO consequence by var id');
ok($var->display_consequence('label') eq 'intron variant', 'display consequence label by var id');
ok( !$var->is_failed(),              "failed status");

my %syns = map {$_ => 1} @{$var->get_all_synonyms()};

ok( $syns{rs57302278}, 'archive rs synonym' );

ok($var->get_all_synonym_sources->[0] eq 'Archive dbSNP', 'synonym source');

my $states = $var->get_all_evidence_values();

ok( $states->[0] eq 'Multiple_observations' &&
    $states->[1] eq 'Frequency' && 
    $states->[2] eq 'HapMap' &&
    $states->[3] eq '1000Genomes', 'evidence status');
ok( $var->add_evidence_value("Cited"), 'add a permitted evidence value' );

my %alleles = map {$_->dbID() => $_} @{$var->get_all_Alleles()};

ok($alleles{228265191}->allele() eq 'T' &&
   $alleles{228265243}->allele() eq 'A' &&
   $alleles{228265243}->population()->name() eq '1000GENOMES:pilot_1_YRI_low_coverage_panel', "allele by id");



# test fetch by name
$var = $va->fetch_by_name('rs142276873');

ok($var->name() eq 'rs142276873', "name by name");
ok($var->dbID() == 30220007,      "id by name" );
ok($var->source_name() eq 'dbSNP',"source by name");
ok($var->ancestral_allele eq 'G', "ancestral allele by name");

# test fetch by subsnp
my $var_ss = $va->fetch_by_subsnp_id('ss11455892');
ok($var_ss->name() eq $var_name, 'fetch by subsnp'); 


# test fetch by name using a synonym

$var = $va->fetch_by_name('rs57302278');

ok($var->name() eq $var_name,   "current name by synonym");
ok($var->dbID() == $var_id,       "current id by synonym");
ok($var->source_name() eq 'dbSNP',"source by synonym");
ok($var->ancestral_allele eq 'A', "ancestral allele by synonym" );



#test ambig_code and var_class - fix core test db & re-install
#ok($var->ambig_code eq 'W',   "ambiguity code by synonym");
ok($var->var_class() eq 'SNP',  "variation class by synonym");

ok($var->get_all_synonym_sources->[0] eq 'Archive dbSNP', "synonym source by synonym");

## check a failed one
$va->db->include_failed_variations(1);
my @fails = ('None of the variant alleles match the reference allele', 'Mapped position is not compatible with reported alleles');
my $fail_desc = join(";", @fails);
my $failed_var = $va->fetch_by_name('rs67521280');
ok($failed_var->name() eq 'rs67521280', "name by name");
ok($failed_var->failed_description() eq $fail_desc,   "fail description"); 
ok(join(";", @{$failed_var->get_all_failed_descriptions()}) eq $fail_desc,   "all fail descriptions");


## Iterators
print "\n## Test - Iterators ##\n";
my %var_list = ( 'rs80359159'  => 20949482,
                 'rs117161559' => 25992950,
                 'rs138574806' => 27085574,
                 'rs183577856' => 40763453
               );
my @var_dbIDs = values(%var_list);
my @var_names = sort{ $var_list{$a} <=> $var_list{$b}} keys(%var_list);

# test fetch Iterator by dbID list
print "\n# Test - fetch_Iterator_by_dbID_list\n";
my $it1 = $va->fetch_Iterator_by_dbID_list(\@var_dbIDs);
ok($var_list{$it1->next()->name}, "iterator by id - 1");
ok($var_list{$it1->next()->name}, "iterator by id - 2");
ok($var_list{$it1->next()->name}, "iterator by id - 3");
ok($var_list{$it1->next()->name}, "iterator by id - 4");

# test fetch Iterator by VariationSet
print "\n# Test - fetch_Iterator_by_VariationSet\n";
my $vs  = $vsa->fetch_by_name('1000 Genomes - AFR');
my $isv = $va->fetch_Iterator_by_VariationSet($vs);
ok($var_list{$isv->next()->name}, "iterator by VariationSet - 1");
ok($var_list{$isv->next()->name}, "iterator by VariationSet - 2");
ok($var_list{$isv->next()->name}, "iterator by VariationSet - 3");
ok($var_list{$isv->next()->name}, "iterator by VariationSet - 4");

my $it2 = $va->fetch_Iterator;
ok($it2 && $it2->isa('Bio::EnsEMBL::Utils::Iterator'), 'is Iterator');

my $it3 = $va->fetch_Iterator_somatic;
ok($it3 && $it3->isa('Bio::EnsEMBL::Utils::Iterator'), 'is Iterator');

print "\nfetch\n";
# fetch_all
my $vars = $va->fetch_all;
ok( scalar (grep { $_->name eq 'rs2255888' } @$vars) == 1, 'fetch_all');

# fetch_all_somatic
$vars = $va->fetch_all_somatic;
ok(scalar @$vars == 1, 'fetch_all_somatic');

# test fetch by stable_id
print "\n# Test - fetch_by_stable_id\n";
my $var2 = $va->fetch_by_stable_id($var_name);
ok($var2->name eq $var_name, "var by stable_id");

# test fetch all by source
print "\n# Test - fetch_all_by_source\n";
my $var3 = $va->fetch_all_by_source('dbSNP');
ok($var3->[0]->source_name eq 'dbSNP', "var by source");

# test fetch all by source type
print "\n# Test - fetch_all_by_source_type\n";
my $var3a = $va->fetch_all_by_source_type('lsdb');
ok($var3a->[0]->name eq 'rs121908760', "var by source type");

# test fetch all by dbID list
print "\n# Test - fetch_all_by_dbID_list\n";
@var_dbIDs = values(%var_list);
my $var4 = $va->fetch_all_by_dbID_list(\@var_dbIDs);
ok($var_list{$var4->[0]->name}, "var by list of dbIDs");

# test fetch all by name list
print "\n# Test - fetch_all_by_name_list\n";
my $var5 = $va->fetch_all_by_name_list(\@var_names);
ok($var_list{$var5->[0]->name}, "var by list of names");

# test get all sources
print "\n# Test - get_all_sources\n";
my $src1 = $va->get_all_sources();
ok($src1->[0] eq 'dbSNP',  "get all sources - 1");
ok($src1->[1] eq 'COSMIC', "get all sources - 2");

# test get flanking sequence
print "\n# Test - get_flanking_sequence\n";
my $fs = $va->get_flanking_sequence($var_id, 10);
ok($fs->[0] eq 'NNNNNNNNNNN', "get the flanking sequence");

# test fetch all by Population
print "\n# Test - fetch_all_by_Population\n";
my $pop = $pa->fetch_by_name('CSHL-HAPMAP:HAPMAP-ASW');
my $var6 = $va->fetch_all_by_Population($pop);
ok( scalar (grep { $_->name eq 'rs2255888' } @$var6) == 1, "var by population");

# test fetch all by VariationSet
print "\n# Test - fetch_all_by_VariationSet\n";
my $var7 = $va->fetch_all_by_VariationSet($vs);
ok($var_list{$var7->[0]->name}, "var by VariationSet - 1");
ok($var_list{$var7->[1]->name}, "var by VariationSet - 2");
ok($var_list{$var7->[2]->name}, "var by VariationSet - 3");

# test load alleles
print "\n# Test - load_alleles\n";
$va->load_alleles(1);
my $var8 = $va->fetch_by_dbID($var_id);
ok($var8->get_all_Alleles, "load alleles");

# store
print "\n# Test - store\n";
$var = $va->fetch_by_dbID($var_id);

delete $var->{$_} for qw(dbID name);
$var->name('test');
$var->add_synonym('dbSNP', 'ss55331');
ok($va->store($var), "store");

$var = $va->fetch_by_name('test');
ok($var && $var->name eq 'test', "fetch stored");
ok($var->get_all_synonyms('dbSNP')->[0] eq 'ss55331', "fetch synonym stored with variant");

# update
print "\n# Test - update\n";
my $upd_name = 'updated_test';
my $upd_var = $va->fetch_by_name('test');
$upd_var->name($upd_name);
$va->update($upd_var);
$var = $va->fetch_by_name($upd_name);
ok($var && $var->name eq $upd_name, "fetch updated");


## attrib handling
$var->update_attributes( {"co-located allele"  => "colo"} );
ok($var->get_all_attributes()->{"co-located allele"}  eq "colo",  "attribute extraction");

$va->store_attributes($var);
my $var_up = $va->fetch_by_dbID($var->dbID);
ok($var_up && $var_up->get_all_attributes()->{"co-located allele"} eq "colo", "fetch updated attribs");

## test synonyms
$var->add_synonym('dbSNP', 'ss55331');
$va->store_synonyms($var);

my $varup = $va->fetch_by_name($upd_name);
ok($varup->get_all_synonyms('dbSNP')->[0] eq 'ss55331', "fetch updated synonym");

## test bad synonym source
$var->add_synonym('turnip', 'ss55331');
throws_ok { $va->store_synonyms($var) } qr/No source found for name turnip/, 'Throw if source not found.';


done_testing();
