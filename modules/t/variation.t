use strict;
use warnings;

BEGIN { $| = 1;
	use Test;
	plan tests => 36;
}


use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Variation::Variation;
use Bio::EnsEMBL::Variation::Allele;

use Bio::EnsEMBL::Test::MultiTestDB;

# test constructor

my $dbID = 123;
my $name = 'rs5432';
my $source = 'dbSNP';
my $synonyms = {'dbSNP' => ['ss355', 'ss556'], 'TSC' => ['12565']};
my $a1 = Bio::EnsEMBL::Variation::Allele->new(-allele => 'A');
my $a2 = Bio::EnsEMBL::Variation::Allele->new(-allele => 'C');
my $alleles = [$a1,$a2];
my $validation_states = ['submitter', 'cluster'];
my $five_prime_seq = 'AAATTAACCATTGGCG';
my $three_prime_seq = 'TTATTTTAAGGCCGGAGTA';
my $ancestral_allele = 'A';
my $moltype = 'Genomic';


my $v = Bio::EnsEMBL::Variation::Variation->new
  (-dbID => 123,
   -name => $name,
   -source => $source,
   -synonyms => $synonyms,
   -alleles => $alleles,
   -validation_states => $validation_states,
   -five_prime_flanking_seq => $five_prime_seq,
   -three_prime_flanking_seq => $three_prime_seq,
   -ancestral_allele => $ancestral_allele,
   -moltype => $moltype);

ok($v->dbID());
ok($v->name() eq $name);
ok($v->source() eq $source);


ok(@{$v->get_all_synonyms()} == 3);
ok($v->get_all_synonyms('TSC')->[0] eq '12565');
ok($v->get_all_Alleles()->[0]->allele() eq 'A');
ok($v->get_all_validation_states()->[0] eq 'cluster' &&
   $v->get_all_validation_states()->[1] eq 'submitter');

ok($v->five_prime_flanking_seq() eq $five_prime_seq);
ok($v->three_prime_flanking_seq() eq $three_prime_seq);

#test amibg_code
ok($v->ambig_code() eq 'M');

#test variation_class
ok($v->var_class() eq 'snp');

##test ancestral_allele
ok($v->ancestral_allele() eq 'A');

##test molecular type
ok($v->moltype() eq 'Genomic');

# test getter/setters

ok(test_getter_setter($v, 'name', 'newname'));
ok(test_getter_setter($v, 'source', 'newsource'));
ok(test_getter_setter($v, 'five_prime_flanking_seq', 'AATTTA'));
ok(test_getter_setter($v, 'three_prime_flanking_seq', 'TTTA'));
ok(test_getter_setter($v,'ancestral_allele','C'));
ok(test_getter_setter($v,'moltype','cDNA'));



# test add_synonym, get_all_synonym_sources and get_all_synonyms

$v->add_synonym('newsource', 'mysyn');
ok($v->get_all_synonyms('newsource')->[0] eq 'mysyn');

my @sources = sort {$a cmp $b} @{$v->get_all_synonym_sources()};

ok($sources[0] eq 'TSC' &&
   $sources[1] eq 'dbSNP' &&
   $sources[2] eq 'newsource');

ok(@{$v->get_all_synonyms()} == 4);



# test add_Allele, get_all_Alleles

my $a3  = Bio::EnsEMBL::Variation::Allele->new('allele' => '-');
$v->add_Allele($a3);
ok($v->get_all_Alleles()->[2] == $a3);


# test add_validation_state
$v->add_validation_state('freq');
# states are always added in same order
ok(join(',', @{$v->get_all_validation_states()}) eq 'cluster,freq,submitter');

# adding the same state twice does nothing
$v->add_validation_state('freq');
# states are always added in same order
ok(join(',', @{$v->get_all_validation_states()}) eq 'cluster,freq,submitter');


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
   -adaptor => $var_adaptor
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

