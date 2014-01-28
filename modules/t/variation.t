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
BEGIN { $| = 1;
	use Test::More;
	plan tests => 23;
}

use FindBin qw($Bin);

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Variation::Variation;
use Bio::EnsEMBL::Variation::Allele;


## adaptor needed as availabilty checked in Allele.pm
my $reg = 'Bio::EnsEMBL::Registry';
$reg->no_version_check(1); ## version not relevant for test db
$reg->load_all("$Bin/test.ensembl.registry");
my $variation_adaptor    = $reg->get_adaptor('homo_sapiens', 'variation', 'variation');


# test constructor

my $dbID = 123;
my $name = 'rs5432';
my $source = 'dbSNP';
my @synonyms = ( ['dbSNP', 'ss355',  1,],
                 ['dbSNP', 'ss556',  1,],
                 ['TSC',   '12565',  1 ] );
my %syonym;
foreach my $synonyms (@synonyms){
    $syonym{$synonyms->[0]}{$synonyms->[1]}{$synonyms->[2]}++;
}

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

ok($v->dbID() eq 123,             "db ID");
ok($v->name() eq $name,           "name");
ok($v->source() eq $source,       "source");
ok($v->is_somatic() eq 0,         "is_somatic");
ok($v->flipped() eq 0,            "flipped");
ok($v->moltype() eq 'Genomic',    "moltype");
ok($v->ancestral_allele() eq 'A', "ancestral_allele");


my $n = scalar @{$v->get_all_synonyms()};

ok(@{$v->get_all_synonyms()} == 3, "count synonym");
ok($v->get_all_synonyms('TSC')->[0] eq '12565', "syonym by source");
ok($v->get_all_Alleles()->[0]->allele() eq 'A', "allele");
ok($v->get_all_validation_states()->[0] eq 'cluster' &&
   $v->get_all_validation_states()->[1] eq 'submitter', "validation");

#test amibg_code
#ok($v->ambig_code() eq 'M', "ambig code");

#test variation_class
ok($v->var_class() eq 'SNP', "class");


# test getter/setters

ok(test_getter_setter($v, 'name', 'newname'), "get/set name");
ok(test_getter_setter($v, 'source', 'newsource'), "get/set source");
ok(test_getter_setter($v, 'source_url', 'http://www.ncbi.nlm.nih.gov/projects/SNP/'), "get/set source url");
ok(test_getter_setter($v, 'ancestral_allele','C'), "get/set ancestral_allele");
ok(test_getter_setter($v, 'moltype','cDNA'), "get/set moltype");
ok(test_getter_setter($v, 'clinical_significance', $clin_sig), "get/set clin_sig");


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


=head test requiring MultiTestDB
my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();

my $vdb = $multi->get_DBAdaptor('variation');
my $core = $multi->get_DBAdaptor('core');
$vdb->dnadb($core);

my $var_adaptor = $vdb->get_VariationAdaptor;
#test get_all_IndividualGenotypes

my $variation_id = 191;

my $variation = Bio::EnsEMBL::Variation::Variation->new(
   -dbID => $variation_id,
   -name => 'rs193',
   -adaptor => $variation_adaptor 
   );

my $igty = $variation->get_all_IndividualGenotypes();
my @igtys = sort {$a->individual->dbID() <=> $b->individual->dbID()}
            @{$variation->get_all_IndividualGenotypes()};
ok(@igtys == 96);
ok($igtys[0]->variation()->name() eq 'rs193');
ok($igtys[0]->allele1() eq 'C');
ok($igtys[0]->allele2() eq 'T');
ok($igtys[0]->individual()->name() eq 'NA17011');

#test get_all_PopulationGenotypes
$variation_id = 2863;

$variation = Bio::EnsEMBL::Variation::Variation->new(
   -dbID => $variation_id,
   -name => 'rs2872',
   -adaptor => $var_adaptor
   );

@igtys = ();

@igtys = sort {$a->dbID() <=> $b->dbID()}
            @{$variation->get_all_PopulationGenotypes()};

ok(@igtys == 12);
ok($igtys[0]->dbID() == 1);
ok($igtys[0]->population()->name() eq 'AFFY:AfAm');
ok($igtys[0]->allele1() eq 'C');
ok($igtys[0]->allele2() eq 'C');
ok($igtys[0]->frequency() == 0.666667);
=cut


done_testing();
