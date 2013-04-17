use strict;
use warnings;

use Test::More;
use FindBin qw($Bin);

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::Population;

my $reg = 'Bio::EnsEMBL::Registry';
$reg->load_all("$Bin/test.ensembl.registry.72");

my $vdba = $reg->get_DBAdaptor('human', 'variation');
my $pa = $vdba->get_PopulationAdaptor;
my $ia = $vdba->get_IndividualAdaptor;
my $va = $vdba->get_VariationAdaptor;
my $vfa = $vdba->get_VariationFeatureAdaptor;
my ($populations, $population, $all);
# fetch_all_1KG_Populations
$populations = $pa->fetch_all_1KG_Populations;
is(scalar @$populations, 31, "Number of 1000 genomes populations");
# fetch_all_HapMap_Populations
$populations = $pa->fetch_all_HapMap_Populations;
is(scalar @$populations, 12, "Number of HapMap populations");
# fetch_all_LD_Populations
$populations = $pa->fetch_all_LD_Populations;
is(scalar @$populations, 26, "Number of LD populations");
# fetch_default_LDPopulation
$population = $pa->fetch_default_LDPopulation;
is($population->name, 'CSHL-HAPMAP:HapMap-CEU', "Name for default LD population");

# fetch_by_dbID
$population = $pa->fetch_by_dbID(145);
is($population->name, 'MPJ6:E-1', "Fetch by dbID 145");
# fetch_all_by_dbID_list
my $list = [145, 400892];
$populations = $pa->fetch_all_by_dbID_list($list);
$all = join(',', map{$_->name} sort {$a->name cmp $b->name} @$populations);
is($all, 'COSMIC:gene:ZZZ3:tumour_site:ovary,MPJ6:E-1', "Fetch by list [145, 400892]");

# fetch_all_by_Individual
# 1000GENOMES:phase_1:HG01625 101495 -> 1000GENOMES:phase_1_ALL,1000GENOMES:phase_1_EUR,1000GENOMES:phase_1_IBS
# CH18564 17630
my @individual_dbids = (101495, 17630);
my $individual = $ia->fetch_by_dbID($individual_dbids[0]);
$populations = $pa->fetch_all_by_Individual($individual);
is(scalar @$populations, 3, "Number of populations for individual HG01625");

# fetch_all_by_Individual_list
# 101495, 17630 -> 1000GENOMES:phase_1_ALL,1000GENOMES:phase_1_EUR,1000GENOMES:phase_1_IBS,CSHL-HAPMAP:HAPMAP-CHB,CSHL-HAPMAP:HapMap-HCB,EGP_SNPS:CHB_GENO_PANEL,ILLUMINA:CHB,PGA-UW-FHCRC:CHB_GENO_PANEL
my $individuals = [];
foreach my $dbid (@individual_dbids) {
   push @$individuals, $ia->fetch_by_dbID($dbid); 
}
$populations = $pa->fetch_all_by_Individual_list($individuals);
is(scalar @$populations, 8, "Number of populations for individuals HG01625 and CH18564");

# fetch_by_name
$population = $pa->fetch_by_name('1000GENOMES:phase_1_IBS');
is($population->name, '1000GENOMES:phase_1_IBS', "Fetch by name 1000GENOMES:phase_1_IBS");

# fetch_all_by_name_search
$populations = $pa->fetch_all_by_name_search('1000GENOMES');
is(scalar @$populations, 31, "Number of populations for fetch by name search = 1000GENOMES");

# fetch_all_by_sub_Population
$population = $pa->fetch_by_name('1000GENOMES:phase_1_IBS');
$populations = $pa->fetch_all_by_sub_Population($population);
$all = join(',', map {$_->name} sort {$a->name cmp $b->name} @$populations);
is($all, '1000GENOMES:phase_1_EUR', "Fetch by sub-population 1000GENOMES:phase_1_IBS");

# fetch_all_by_super_Population
$population = $pa->fetch_by_name('1000GENOMES:phase_1_EUR');
# 1000GENOMES:phase_1_CEU,1000GENOMES:phase_1_FIN,1000GENOMES:phase_1_GBR,1000GENOMES:phase_1_IBS,1000GENOMES:phase_1_TSI
$populations = $pa->fetch_all_by_super_Population($population);
is(scalar @$populations, 5, "Fetch by super-population 1000GENOMES:phase_1_EUR");

# fetch_population_by_synonym
$populations = $pa->fetch_population_by_synonym(627);
$all = join(',', map {$_->name} sort {$a->name cmp $b->name} @$populations);
is($all, 'MPJ6:E-1', "Fetch by synonym 627");

# fetch_synonyms
my $synonyms = $pa->fetch_synonyms(145);
is(join(',', @$synonyms), 627, "Fetch synonyms for 145");

# fetch_tag_Population
#my $v = $va->fetch_by_name('rs205621');
#my $vfs = $vfa->fetch_all_by_Variation($v);
#$populations = $pa->fetch_tag_Population($vfs->[0]);

# fetch_tagged_Population

# get_sample_id_for_population_names

done_testing();
