use strict;
use warnings;

BEGIN { $| = 1;
	use Test;
	plan tests => 8;
}


use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Variation::Allele;
use Bio::EnsEMBL::Variation::Population;




# test constructor

my $dbID = 1;
my $allele = 'A';
my $frequency = 0.86;
my $p = Bio::EnsEMBL::Variation::Population->new();

my $al = Bio::EnsEMBL::Variation::Allele->new
  (-dbID => $dbID,
   -allele => $allele,
   -frequency => $frequency,
   -population => $p);

ok($al->dbID() == $dbID);
ok($al->frequency() == $frequency);
ok($al->population() == $p);
ok($al->allele() eq $allele);



# test getter/setters

my $p2 = Bio::EnsEMBL::Variation::Population->new();

ok(test_getter_setter($al, 'dbID', 123));
ok(test_getter_setter($al, 'allele', 'T'));
ok(test_getter_setter($al, 'frequency', 0.86));
ok(test_getter_setter($al, 'population', $p2));
