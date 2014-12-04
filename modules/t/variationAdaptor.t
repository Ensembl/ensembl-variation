# Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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


use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Test::MultiTestDB;


our $verbose = 0;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');

my $vdb = $multi->get_DBAdaptor('variation');

my $va = $vdb->get_VariationAdaptor();

ok($va && $va->isa('Bio::EnsEMBL::Variation::DBSQL::VariationAdaptor'), "isa var adaptor");

# test fetch by dbID

my $var = $va->fetch_by_dbID(4770800);

ok($var->name() eq 'rs7569578',      'variation name by var id');
ok($var->stable_id() eq 'rs7569578',      'stable id by var id');
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
ok($var->display_consequence() eq 'Intron variant', 'display status by var id');
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
ok($var_ss->name() eq 'rs7569578', 'fetch by subsnp'); 


# test fetch by name using a synonym

$var = $va->fetch_by_name('rs57302278');

ok($var->name() eq 'rs7569578',   "current name by synonym");
ok($var->dbID() == 4770800,       "current id by synonym");
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

done_testing();
