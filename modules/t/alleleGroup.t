use strict;
use warnings;

BEGIN { $| = 1;
	use Test;
	plan tests => 14;
}

use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Variation::AlleleGroup;
use Bio::EnsEMBL::Variation::VariationGroup;
use Bio::EnsEMBL::Variation::Population;
use Bio::EnsEMBL::Variation::Variation;


# test constructor

my $dbID = 4;
my $name = 'test allele grp';
my $pop = Bio::EnsEMBL::Variation::Population->new(-name => 'test');

my $src = 'dbSNP';
my $freq = 0.85;

my $vg = Bio::EnsEMBL::Variation::VariationGroup->new(-name => 'test vg');

my $ag = Bio::EnsEMBL::Variation::AlleleGroup->new
  (-dbID => $dbID,
   -name => $name,
   -population  => $pop,
   -source => $src,
   -variation_group => $vg,
   -frequency => $freq);


ok($ag->dbID() == $dbID);
ok($ag->name() eq $name);
ok($ag->population() == $pop);
ok($ag->source() eq $src);
ok($ag->variation_group() == $vg);
ok($ag->frequency() == $freq);


# test add_Variation , get_all_Variations and get_all_alleles
my $var1 = Bio::EnsEMBL::Variation::Variation->new(-name => 'var1');
my $var2 = Bio::EnsEMBL::Variation::Variation->new(-name => 'var2');

$ag->add_Variation($var1, 'A');
$ag->add_Variation($var2, 'T');

my $vars = $ag->get_all_Variations();
my $alleles = $ag->get_all_alleles();

ok($vars->[0]->name() eq 'var1');
ok($vars->[1]->name() eq 'var2');

ok($alleles->[0] eq 'A');
ok($alleles->[1] eq 'T');



# test getter/setters
ok(test_getter_setter($ag, 'name', 'new name'));
my $new_pop = Bio::EnsEMBL::Variation::Population->new(-name => 'new pop');
ok(test_getter_setter($ag, 'population', $new_pop));
my $new_vg = Bio::EnsEMBL::Variation::VariationGroup->new(-name => 'new vg');
ok(test_getter_setter($ag, 'variation_group', $new_vg));
ok(test_getter_setter($ag, 'frequency', 0.21));





1;
