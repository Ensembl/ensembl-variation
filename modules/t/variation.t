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
use Bio::EnsEMBL::Variation::Population;
use Bio::EnsEMBL::Variation::PopulationGenotype;


my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');

my $vdb = $multi->get_DBAdaptor('variation');
my $core = $multi->get_DBAdaptor('core');

my $variation_adaptor = $vdb->get_VariationAdaptor;

# test constructor

## need source object 

my $source_name           = 'dbSNP';
my $source_version        = 138;
my $source_description    = 'Variants (including SNPs and indels) imported from dbSNP (mapped to GRCh38)';
my $source_type           = 'archive';
my $source_somatic_status = 'mixed';

my $source = Bio::EnsEMBL::Variation::Source->new
  (-name           => $source_name,
   -version        => $source_version,
   -description    => $source_description,
   -type           => $source_type,
   -somatic_status => $source_somatic_status
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
my $ancestral_allele = 'A';
my $clin_sig = 'untested';


my $v = Bio::EnsEMBL::Variation::Variation->new
  (-dbID              => 123,
   -name              => $name,
   -source            => $source,
   -synonyms          => \%syonym,
   -alleles           => $alleles,
   -adaptor           => $variation_adaptor,
   -ancestral_allele  => $ancestral_allele,
   -is_somatic        => 0,
   -flipped           => 0,
  );

ok($v->dbID() eq 123,             "db ID");
ok($v->name() eq $name,           "name");
ok($v->is_somatic() eq 0,         "is_somatic");
ok($v->flipped() eq 0,            "flipped");
ok($v->ancestral_allele() eq 'A', "ancestral_allele");
# source
ok($v->source->name() eq $source_name,                  'source' );
ok($v->source_name eq $source_name,                     'source_name');
ok($v->source_description eq $source_description,       'source_description');
ok($v->source_version eq $source_version,               'source_version');
ok($v->source_type eq $source_type,                     'source_type');
ok($v->source_somatic_status eq $source_somatic_status, 'source_somatic_status');
ok($v->has_somatic_source == 0,                         'has_somatic_source');
ok($v->has_failed_alleles == 0,                         'has_failed_alleles');

# test source object
ok($v->source($source), 'source (using argument)');

ok($v->five_prime_flanking_seq('ATGCCC') eq 'ATGCCC',    'get/set flank');
ok($v->three_prime_flanking_seq('CTTGTA') eq 'CTTGTA',    'get/set flank');

my $n = scalar @{$v->get_all_synonyms()};

ok(@{$v->get_all_synonyms()} == 3, "count synonym");
ok($v->get_all_synonyms('TSC')->[0] eq '12565', "synonym by source");
ok($v->get_all_Alleles()->[0]->allele() eq 'A', "allele");

#test amibg_code
ok($v->ambig_code() eq 'M', "ambig code");

#test variation_class
ok($v->var_class() eq 'SNP', "class");

# test variation_class without alleles
my $v2 = $variation_adaptor->fetch_by_dbID(112643842);
delete($v2->{$_}) for qw(class_display_term class_SO_term);
is($v2->var_class(), 'SNP', 'class via variation feature');


# test getter/setters
ok(test_getter_setter($v, 'name', 'newname'), "get/set name");
ok(test_getter_setter($v, 'ancestral_allele','C'), "get/set ancestral_allele");


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


#test get_all_SampleGenotypes

my $variation = $variation_adaptor->fetch_by_dbID(1748253);

my $sgty = $variation->get_all_SampleGenotypes();
my @sgtys = sort {$a->sample->dbID() <=> $b->sample->dbID()}
            @{$variation->get_all_SampleGenotypes()};

ok(@sgtys == 2,                                   "sample geno to count" );
ok($sgtys[0]->variation()->name() eq 'rs2299222', "sample geno to var name");
ok($sgtys[0]->allele1() eq 'T',                   "sample geno to allele1");
ok($sgtys[0]->allele2() eq 'T',                   "sample geno to allele2");
ok($sgtys[0]->sample()->name() eq 'NA12891',  "sample geno to dna name");



# test get_all_PopulationGenotypes
my $variation_p =  $variation_adaptor->fetch_by_name('rs17081232');
my $pop = $variation_adaptor->db->get_PopulationAdaptor()->fetch_by_name('PERLEGEN:AFD_EUR_PANEL');

my @pgtys = sort {$a->dbID() <=> $b->dbID()} @{$variation_p->get_all_PopulationGenotypes()};

ok(@pgtys == 14,                                  "pop geno to count" );
ok($pgtys[0]->dbID() == 50885783,                 "pop geno to dbID");
ok($pgtys[0]->allele1() eq 'A',                   "pop geno to allele1");
ok($pgtys[0]->allele2() eq 'G',                   "pop geno to allele1");
ok($pgtys[0]->frequency() == 0.125,               "pop geno to freq");
ok($pgtys[0]->population()->name() eq $pop->name,"pop geno to pop name");

# test get_all_Alleles($pop)
my @p_alleles = @{$variation_p->get_all_Alleles($pop) || []};
my $total_freq = 0;
$total_freq += $_->frequency for @p_alleles;

ok(@p_alleles == 2,                                             "get alleles by pop - count");
ok($p_alleles[0]->population->name eq 'PERLEGEN:AFD_EUR_PANEL', "get alleles by pop - name");
ok(join(",", sort map {$_->allele()} @p_alleles) eq 'A,G',      "get alleles by pop - alleles");
ok($total_freq == 1,                                            "get alleles by pop - total frequency");

# test get_all_PopulationGenotypes($pop)
my @p_pgtys = @{$variation_p->get_all_PopulationGenotypes($pop) || []};
$total_freq = 0;
$total_freq += $_->frequency for @p_pgtys;

ok(@p_alleles == 2,                                           "get pop geno by pop - count");
ok($p_pgtys[0]->population->name eq $pop->name, "get pop geno by pop - name");
ok($total_freq == 1,                                          "get pop geno by pop - total frequency");

# derived allele
my $da = $variation_p->derived_allele($pop);
is($da, 'A', "derived allele");

# derived allele freq
my $df = $variation_p->derived_allele_frequency($pop);
is($df, '0.0625', "derived allele frequency");


# test add_PopulationGenotype

my $new_pop = Bio::EnsEMBL::Variation::Population->new(-name => 'New population');
my $new_pg  = Bio::EnsEMBL::Variation::PopulationGenotype->new(
  -genotype      => ['A','A'],
  -frequency     => 0.4,
  -population    => $new_pop,
  -variation     => $variation_p,
  -_variation_id => $variation_p->dbID
);
ok($variation_p->add_PopulationGenotype($new_pg), 'add_PopulationGenotype');


# get all Genes
my @genes = @{$variation_adaptor->fetch_by_name('tmp__')->get_all_Genes};
ok(scalar @genes == 1, "get_all_Genes count");
is($genes[0]->stable_id, 'ENSG00000236003', "get_all_Genes stable_id");


# get all clinical significance states
my $var_cs = $variation_adaptor->fetch_by_name('rs121908760');
ok($var_cs->get_all_clinical_significance_states()->[0] eq 'pathogenic', 'get_all_clinical_significance_states');


# get all PhenotypeFeatures
my $pfs = $variation->get_all_PhenotypeFeatures();
ok($pfs->[0]->phenotype->description eq 'ACHONDROPLASIA', 'get_all_PhenotypeFeatures');

# get all Publications
my $var_cit = $variation_adaptor->fetch_by_name('rs7698608');
my $pubs = $var_cit->get_all_Publications();
ok($pubs->[0]->pmid == 22779046, 'get_all_Publications');

## Variation Feature ##

my $var_name = 'rs142276873';
my $var = $variation_adaptor->fetch_by_name($var_name);

# ambiguity code
ok($var->ambig_code() eq 'R', "ambig code, using variation_feature");





# test get_all_VariationFeatures


# test get_VariationFeature_by_dbID
my $vf2 = $var->get_VariationFeature_by_dbID(33303674);
ok($vf2->variation_name eq $var_name, 'get_VariationFeature_by_dbID');


my $hash = {dbID              => 123,
            name              => $name,
            source            => $source,
            synonyms          => \%syonym,
            alleles           => $alleles,
            adaptor           => $variation_adaptor
           };

my $vfast =  Bio::EnsEMBL::Variation::Variation->new_fast($hash);

ok($vfast->dbID() eq 123,             "fast - db ID");
ok($vfast->name() eq $name,           "fast - name");
ok(! $vfast->ancestral_allele()  , "no derived");


done_testing();
