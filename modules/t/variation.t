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
use Data::Dumper;
use Test::More;
use Bio::EnsEMBL::Test::MultiTestDB;

use FindBin qw($Bin);

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Variation::Variation;
use Bio::EnsEMBL::Variation::Allele;
use Bio::EnsEMBL::Variation::Source;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');

my $vdb = $multi->get_DBAdaptor('variation');
my $core = $multi->get_DBAdaptor('core');

my $variation_adaptor = $vdb->get_VariationAdaptor;


# test constructor

## need source object 

my $source_name           = 'dbSNP';
my $source_version        = 138;
my $source_description    = 'Variants (including SNPs and indels) imported from dbSNP (mapped to GRCh38)';

my $source = Bio::EnsEMBL::Variation::Source->new
  (-name           => $source_name,
   -version        => $source_version,
   -description    => $source_description
);


my $dbID = 123;
my $name = 'rs5432';
my @synonyms = ( ['dbSNP', 'ss355',  1,],
                 ['dbSNP', 'ss556',  1,],
                 ['TSC',   '12565',  1 ] );
my %syonym;
foreach my $synonyms (@synonyms){
    $syonym{$synonyms->[0]}{$synonyms->[1]}{$synonyms->[2]}++;
}

## need allele objects
my $a1 = Bio::EnsEMBL::Variation::Allele->new(-allele => 'A', -adaptor => $variation_adaptor);
my $a2 = Bio::EnsEMBL::Variation::Allele->new(-allele => 'C', -adaptor => $variation_adaptor);
my $alleles = [$a1,$a2];
my $validation_states = ['submitter', 'cluster'];
my $ancestral_allele = 'A';
my $moltype = 'Genomic';
my $clin_sig = 'untested';


my $v = Bio::EnsEMBL::Variation::Variation->new
  (-dbID              => 123,
   -name              => $name,
   -source            => $source,
   -synonyms          => \%syonym,
   -alleles           => $alleles,
   -adaptor           => $variation_adaptor,
   -validation_states => $validation_states,
   -ancestral_allele  => $ancestral_allele,
   -moltype           => $moltype,
   -is_somatic        => 0,
   -flipped           => 0,
  );

ok($v->dbID() eq 123,                "db ID");
ok($v->name() eq $name,              "name");
ok($v->source_name() eq $source_name,"source");
ok($v->is_somatic() eq 0,            "is_somatic");
ok($v->flipped() eq 0,               "flipped");
ok($v->moltype() eq 'Genomic',       "moltype");
ok($v->ancestral_allele() eq 'A',    "ancestral_allele");


my $n = scalar @{$v->get_all_synonyms()};

ok(@{$v->get_all_synonyms()} == 3, "count synonym");
ok($v->get_all_synonyms('TSC')->[0] eq '12565', "syonym by source");
ok($v->get_all_Alleles()->[0]->allele() eq 'A', "allele");
ok($v->get_all_validation_states()->[0] eq 'cluster' &&
   $v->get_all_validation_states()->[1] eq 'submitter', "validation");

#test amibg_code
ok($v->ambig_code() eq 'M', "ambig code");

#test variation_class
ok($v->var_class() eq 'SNP', "class");


# test getter/setters

ok(test_getter_setter($v, 'name', 'newname'), "get/set name");
ok(test_getter_setter($v, 'ancestral_allele','C'), "get/set ancestral_allele");
ok(test_getter_setter($v, 'moltype','cDNA'), "get/set moltype");



# test add_synonym and get_all_synonym_sources

$v->add_synonym('newsource', 'mysyn');
my @sources = sort {$a cmp $b} @{$v->get_all_synonym_sources()};

ok($sources[0] eq 'TSC' &&
   $sources[1] eq 'dbSNP' &&
   $sources[2] eq 'newsource', "synonyms");


my @dbsnp_syns = @{$v->get_all_synonyms('dbSNP')};
ok(@{$v->get_all_synonyms()} == 4, "count synonyms");



# test add_Allele, get_all_Alleles

my $a3  = Bio::EnsEMBL::Variation::Allele->new( -allele => '-', -adaptor => $variation_adaptor);
$v->add_Allele($a3);
ok($v->get_all_Alleles()->[2] == $a3, "adding allele");


# test add_validation_state
$v->add_validation_state('freq');
# states are always added in same order
ok(join(',', @{$v->get_all_validation_states()}) eq 'cluster,freq,submitter', "valiation states 1");

# adding the same state twice does nothing
$v->add_validation_state('freq');
# states are always added in same order
ok(join(',', @{$v->get_all_validation_states()}) eq 'cluster,freq,submitter', "valiation states 2");




#test get_all_IndividualGenotypes

my $variation = $variation_adaptor->fetch_by_dbID(1748253);

my $igty = $variation->get_all_IndividualGenotypes();
my @igtys = sort {$a->individual->dbID() <=> $b->individual->dbID()}
            @{$variation->get_all_IndividualGenotypes()};

ok(@igtys == 2,                                   "ind geno to count" );
ok($igtys[0]->variation()->name() eq 'rs2299222', "ind geno to var name");
ok($igtys[0]->allele1() eq 'T',                   "ind geno to allele1");
ok($igtys[0]->allele2() eq 'T',                   "ind geno to allele2");
ok($igtys[0]->individual()->name() eq 'NA12891',  "ind geno to dna name");



#test get_all_PopulationGenotypes

my $variation_p =  $variation_adaptor->fetch_by_name('rs17081232');

my @pgtys = sort {$a->dbID() <=> $b->dbID()}
            @{$variation_p->get_all_PopulationGenotypes()};

ok(@pgtys == 14,                                  "pop geno to count" );
ok($pgtys[0]->dbID() == 50885783,                 "pop geno to dbID");
ok($pgtys[0]->allele1() eq 'A',                   "pop geno to allele1");
ok($pgtys[0]->allele2() eq 'G',                   "pop geno to allele1");
ok($pgtys[0]->frequency() == 0.125,               "pop geno to freq");
ok($pgtys[0]->population()->name() eq 'PERLEGEN:AFD_EUR_PANEL',"pop geno to pop name");


done_testing();
