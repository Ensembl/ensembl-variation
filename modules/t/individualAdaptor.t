use strict;
use warnings;

use Test::More;
use FindBin qw($Bin);

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::Individual;

my $reg = 'Bio::EnsEMBL::Registry';
$reg->load_all("$Bin/test.ensembl.registry.72");
my $vdba = $reg->get_DBAdaptor('human', 'variation');

my $pa = $vdba->get_PopulationAdaptor;
my $ia = $vdba->get_IndividualAdaptor;

my ($tests, $individuals, $individual, $all);

# Test IndividualAdaptor
# fetch_all_by_Population 
$tests = [
    { population => '1000GENOMES:phase_1_AFR', size => 246,},
    { population => '1000GENOMES:phase_1_TSI', size => 98,},
    { population => 'CSHL-HAPMAP:HAPMAP-MEX', size => 90,},
];

foreach my $test (@$tests) {
    my $population = $pa->fetch_by_name($test->{population});
    my $individuals = $ia->fetch_all_by_Population($population);
    my $size = scalar @$individuals;
    is($size, $test->{size}, "Individual count for $test->{population}");
}

# fetch_all_by_name
$tests = [
    { name => 'Venter', count => 1,},
    { name => '1000GENOMES', count => 0,},
    { name => 'NA12335', count => 1,},
];

foreach my $test (@$tests) {
    my $individuals = $ia->fetch_all_by_name($test->{name});
    my $count = scalar @$individuals;
    is($count, $test->{count}, "Number of returned individuals for $test->{name}");
}

# fetch_by_dbID
$individual = $ia->fetch_by_dbID(12987);
is($individual->name, 'CEPH104.03', 'Fetch by dbID 12987');

# fetch_all_by_dbID_list -- needed by web team..
my $list = [12987];
$individuals = $ia->fetch_all_by_dbID_list($list);
foreach my $individual (@$individuals) {
    is($individual->name, 'CEPH104.03', "Fetch by dbID list [12987]");
}

# fetch_all_by_parent_Individual
my $parent = $ia->fetch_by_dbID(13949);
is($parent->name, 'CEPH104.02', "Parent name is CEPH104.02");
is($parent->description, 'caucasian.took in ind4421.', "Description for CEPH104.02: caucasian.took in ind4421.");

my $populations = $parent->get_all_Populations();
$all = join(',', map {$_->name} sort {$a->name cmp $b->name} @$populations);
is($all,'AFFY:CEPH,TSC-CSHL:CEL_caucasian_CEPH,TSC-CSHL:MOT_caucasian_CEPH', "Populations membership for CEPH104.02");

my $children = $ia->fetch_all_by_parent_Individual($parent);
$all = join(',', map {$_->name} sort {$a->name cmp $b->name} @$children);
is($all, '904_NA11038,CEPH104.03,CEPH104.04,CEPH104.05,CEPH104.06,CEPH104.07,CEPH104.08,CEPH104.09,CEPH104.10,CEPH104.11,CEPH104.12', "All children for CEPH104.02");

# fetch_individual_by_synonym 
my $synonyms = $ia->fetch_synonyms(12991);
$all = join(',', map {$_} @$synonyms);
is($all, 5, "Synonyms for Individual with dbID 12991"); 
my $synonym = '5';
$individual = $ia->fetch_individual_by_synonym($synonym);
$all = join(',', map {$_->name} @$individuals);
is($all, 'CEPH104.03', "Individuals for synonym 5");

# fetch_all_strains
my $strains = $ia->fetch_all_strains();
is(scalar @$strains, 0, "Number of strains");

# fetch_all_strains_with_coverage:
# 1000GENOMES:pilot_2_trio:NA12878,1000GENOMES:pilot_2_trio:NA12891,1000GENOMES:pilot_2_trio:NA12892,1000GENOMES:pilot_2_trio:NA19238,1000GENOMES:pilot_2_trio:NA19239,1000GENOMES:pilot_2_trio:NA19240,AK1,Anonymous Irish Male,Henry Louis Gates Jr,Henry Louis Gates Sr,Marjolein Kriek,Misha Angrist,Palaeo-Eskimo Saqqaq individual,Palaeo-Eskimo Saqqaq individual HC,Rosalynn Gill,SJK,Stephen Quake,VENTER,WATSON,YH
$strains = $ia->fetch_all_strains_with_coverage();
is(scalar @$strains, 20, "Number of strains with coverage");

# get_default_strains: VENTER,WATSON
$strains = $ia->get_default_strains;
is(scalar @$strains, 2, "Number of default strains");

# get_display_strains
$strains = $ia->get_display_strains;
is(scalar @$strains, 21, "Number of display strains");

# get_reference_strain_name

done_testing();
