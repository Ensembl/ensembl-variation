# Copyright 2013 Ensembl
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

BEGIN { $| = 1;
	use Test;
	plan tests => 36;
}


use Bio::EnsEMBL::Test::TestUtils;


use Bio::EnsEMBL::Test::MultiTestDB;


our $verbose = 0;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();

my $vdb = $multi->get_DBAdaptor('variation');

my $va = $vdb->get_VariationAdaptor();

ok($va && $va->isa('Bio::EnsEMBL::Variation::DBSQL::VariationAdaptor'));

# test fetch by dbID

my $var = $va->fetch_by_dbID(8);

ok($var->name() eq 'rs10');
ok($var->source eq 'dbSNP');
ok($va->get_source_version('dbSNP') == 122);
ok($var->ancestral_allele() eq 'A');
ok($var->moltype() eq 'Genomic');

my %syns = map {$_ => 1} @{$var->get_all_synonyms()};

ok($syns{ss17192627}
   && $syns{ss17950970}
   && $syns{ss20325024}
   && $syns{ss4917294}
   && $syns{ss9});

ok($var->get_all_synonym_sources->[0] eq 'dbSNP');

my $states = $var->get_all_validation_states();

ok($states->[0] eq 'cluster' && $states->[1] eq 'freq' &&
   $states->[2] eq 'submitter' && $states->[3] eq 'hapmap');

my %alleles = map {$_->dbID() => $_} @{$var->get_all_Alleles()};

ok($alleles{603}->allele() eq 'A' &&
   $alleles{604}->allele() eq 'C' &&
   $alleles{605}->allele() eq 'A' &&
   $alleles{606}->allele() eq 'C' &&
   $alleles{606}->population()->name() eq 'EGP_SNPS:PDR90');



# test fetch by name
$var = $va->fetch_by_name('rs3');

ok($var->name() eq 'rs3');
ok($var->dbID() == 1);
ok($var->source eq 'dbSNP');
ok($var->ancestral_allele eq 'C');
ok($var->moltype eq 'cDNA');

my @syns = @{$var->get_all_synonyms()};
ok(@syns == 1 && $syns[0] eq 'ss2');

ok($var->get_all_synonym_sources->[0] eq 'dbSNP');

$states = $var->get_all_validation_states();

ok($states->[0] eq 'freq' && $states->[1] eq 'submitter');

my %a = map {$_->dbID() => $_} @{$var->get_all_Alleles()};

ok($a{1}->allele() eq 'C' && $a{1}->frequency() == 0.96 &&
   $a{2}->allele() eq 'T' && $a{2}->frequency() == 0.04 &&
   $a{3}->allele() eq 'C' && $a{3}->frequency() == 0.89 &&
   $a{4}->allele() eq 'T' && $a{4}->frequency() == 0.11 &&
   $a{5}->allele() eq 'C' && $a{5}->frequency() == 0.95 &&
   $a{6}->allele() eq 'T' && $a{6}->frequency() == 0.05 &&
   $a{7}->allele() eq 'C' && $a{7}->frequency() == 0.93 &&
   $a{8}->allele() eq 'T' && $a{8}->frequency() == 0.07);


ok($a{1}->population->name() eq 'KWOK:C');
ok($a{8}->population->name() eq 'KWOK:S');





# test fetch by name using a synonym

$var = $va->fetch_by_name('ss2');

ok($var->name() eq 'rs3');
ok($var->dbID() == 1);
ok($var->source eq 'dbSNP');
ok($var->ancestral_allele eq 'C');
ok($var->moltype eq 'cDNA');

@syns = @{$var->get_all_synonyms()};
ok(@syns == 1 && $syns[0] eq 'ss2');

ok($var->get_all_synonym_sources->[0] eq 'dbSNP');

$states = $var->get_all_validation_states();

ok($states->[0] eq 'freq' && $states->[1] eq 'submitter');

%a = map {$_->dbID() => $_} @{$var->get_all_Alleles()};

ok($a{1}->allele() eq 'C' && $a{1}->frequency() == 0.96 &&
   $a{2}->allele() eq 'T' && $a{2}->frequency() == 0.04 &&
   $a{3}->allele() eq 'C' && $a{3}->frequency() == 0.89 &&
   $a{4}->allele() eq 'T' && $a{4}->frequency() == 0.11 &&
   $a{5}->allele() eq 'C' && $a{5}->frequency() == 0.95 &&
   $a{6}->allele() eq 'T' && $a{6}->frequency() == 0.05 &&
   $a{7}->allele() eq 'C' && $a{7}->frequency() == 0.93 &&
   $a{8}->allele() eq 'T' && $a{8}->frequency() == 0.07);


ok($a{1}->population->name() eq 'KWOK:C');
ok($a{8}->population->name() eq 'KWOK:S');




# test fetch_by_dbID_list
my $list = [1, 2, 3, 4, 10];

my @vars =  sort {$a->dbID <=> $b->dbID} @{$va->fetch_all_by_dbID_list($list)};

ok(@vars == 5);

ok($vars[0]->dbID() == 1 && $vars[4]->dbID() == 10);



#test ambig_code and var_class
ok($var->ambig_code eq 'Y');

ok($var->var_class() eq 'snp');
