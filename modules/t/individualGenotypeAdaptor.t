use strict;
use warnings;

#BEGIN { $| = 1;
#	use Test;
#	plan tests => 8;
#}

#use Bio::EnsEMBL::Test::TestUtils;
#use Bio::EnsEMBL::Test::MultiTestDB;
#my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
#my $vdb = $multi->get_DBAdaptor('variation');

#my $igty_adaptor = $vdb->get_IndividualGenotypeAdaptor();
#ok($igty_adaptor->isa('Bio::EnsEMBL::Variation::DBSQL::CompressedGenotypeAdaptor'));

use Test::More;
use FindBin qw($Bin);

use Bio::EnsEMBL::Registry;

my $reg = 'Bio::EnsEMBL::Registry';
$reg->load_all("$Bin/test.ensembl.registry.72");

my $cdba = $reg->get_DBAdaptor('human', 'core');
my $vdba = $reg->get_DBAdaptor('human', 'variation');

my $va = $vdba->get_VariationAdaptor();
my $pa = $vdba->get_PopulationAdaptor();
my $ia = $vdba->get_IndividualAdaptor();
my $igta  = $vdba->get_IndividualGenotypeAdaptor();
my $igtfa = $vdba->get_IndividualGenotypeFeatureAdaptor(); 
my $variation_name = 'rs144235347';
my $variation = $va->fetch_by_name($variation_name);
is($variation->name, 'rs144235347', 'variation name is rs144235347');
my $population_name = '1000GENOMES:phase_1_IBS';
my $population = $pa->fetch_by_name($population_name);
is($population->name, '1000GENOMES:phase_1_IBS', 'population name is 1000GENOMES:phase_1_IBS');
my $individuals = $population->get_all_Individuals;
is(scalar @$individuals, 14, 'Number of individuals in 1000GENOMES:phase_1_IBS');

my $sa = $cdba->get_SliceAdaptor;
my $slice = $sa->fetch_by_region('chromosome', 13, 54950000, 55000000);

#my $igts = $igta->fetch_all_by_Slice($slice);
#print 'fetch IGTs by slice: ', scalar @$igts, "\n";
my $igts = $igta->fetch_all_by_Slice($slice, $population);
is(scalar @$igts, 8974, 'Fetch IGTs by slice and population');
$igts = $igta->fetch_all_by_Variation($variation, $population);
is(scalar @$igts, 1092, 'Fetch IGTs by variation and population');
my $igtfs = $igtfa->fetch_all_by_Variation($variation); 
is(scalar @$igtfs, 1092, 'Fetch IGTFs by variation');

$individuals = $ia->fetch_all_by_name('1000GENOMES:phase_1:HG00151'); 
#print 'fetch individual by name 1000GENOMES:phase_1:HG00151: ', scalar @$individuals, "\n";
$igtfs = $igtfa->fetch_all_by_Variation($variation, $individuals->[0]);
is(@$igtfs, 1, 'Fetch IGTFs by variation and individual');

my $adaptor = $igtfs->[0]->adaptor;
ok(defined $adaptor && $adaptor->isa('Bio::EnsEMBL::Variation::DBSQL::IndividualGenotypeFeatureAdaptor'), 'Adaptor is attached');

is($igtfs->[0]->individual->name, '1000GENOMES:phase_1:HG00151', 'Individual name');
is($igtfs->[0]->seq_region_start, 3978985, 'seq_start');
is($igtfs->[0]->seq_region_end, 3978985, 'seq_end');
is($igtfs->[0]->seq_region_strand, 1, 'seq_strand');
is($igtfs->[0]->seqname, 'chromosome:GRCh37:10:1:135534747:1', 'seq_name');
is(join('|', @{$igtfs->[0]->genotype}), 'C|C', 'genotype');
is($igtfs->[0]->genotype_string, 'C|C', 'genotype_string');
is($igtfs->[0]->variation->name, 'rs144235347', 'variation_name');

#print 'subsnp:          ', $igtfs->[0]->subsnp, "\n";
#print 'subsnp_handle:   ', $igtfs->[0]->subsnp_handle, "\n";

# fetch_all_by_Slice
# fetch_all_by_Slice_Individual
# fetch_all_by_Slice_Individual_List
# fetch_all_by_Slice_Population
# fetch_all_by_Variation
# fetch_all_by_Variation_Individual
# fetch_all_by_Variation_Population
# fetch_all_by_Variation_IndividualList

done_testing();


# # test fetch_all_by_individual
# my $ind_adaptor = $vdb->get_IndividualAdaptor();
# my $ind = shift @{$ind_adaptor->fetch_individual_by_synonym(1208)};

# my @igtys = sort {$a->variation->dbID() <=> $b->variation->dbID()}
#             @{$igty_adaptor->fetch_all_by_Individual($ind)};

# ok(@igtys == 17);
# ok($igtys[0]->variation()->name() eq 'rs193');
# ok($igtys[0]->allele1() eq 'C');
# ok($igtys[0]->allele2() eq 'T');
# ok($igtys[0]->individual()->dbID() == 1776);


# test fetch_all_by_Variation
#my $variation_adaptor = $vdb->get_VariationAdaptor();
#my $variation = $variation_adaptor->fetch_by_dbID(191);

#my @igtys = ();

#@igtys = sort {$a->individual->dbID() <=> $b->individual->dbID()} @{$igty_adaptor->fetch_all_by_Variation($variation)};

#ok(@igtys == 50);
#ok($igtys[0]->variation()->name() eq 'rs193');
#ok($igtys[0]->allele1() eq 'C');
#ok($igtys[0]->allele2() eq 'T');
#ok($igtys[0]->individual()->name() eq 'NA17011');

#test get_all_Populations
#my @pops = sort {$a->dbID() <=> $b->dbID()} @{$igtys[0]->individual->get_all_Populations()};
#ok(@pops == 3);
#ok($pops[0]->name eq 'TSC-CSHL:CEL_asian');
