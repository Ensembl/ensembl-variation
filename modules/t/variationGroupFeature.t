use strict;
use warnings;

BEGIN { $| = 1;
	use Test;
	plan tests => 8;
}


use Bio::EnsEMBL::Test::TestUtils;

use Bio::EnsEMBL::Variation::VariationGroup;
use Bio::EnsEMBL::Variation::VariationGroupFeature;


# test constructor

my $vg = Bio::EnsEMBL::Variation::VariationGroup->new(-name => 'PERLGEN:00007',
                                                      -source => 'dbSNP');


my $start = 100;
my $end   = 100;
my $strand = 1;
my $vgname = $vg->name();


my $vgf = Bio::EnsEMBL::Variation::VariationGroupFeature->new
  (-start => 100,
   -end   => 100,
   -strand => 1,
   -variation_group_name => $vgname,
   -variation_group => $vg);

ok($vgf->start()   == $start);
ok($vgf->end()     == $end);
ok($vgf->strand()  == $strand);
ok($vgf->variation_group_name() eq $vgname);
ok($vgf->variation_group() == $vg);


# test getter/setters

my $vg2 = Bio::EnsEMBL::Variation::VariationGroup->new
  (-name => 'PERLGEN:00001245',
   -source => 'dbSNP');

ok(test_getter_setter($vgf, 'variation_group', $vg2));
ok(test_getter_setter($vgf, 'variation_group_name', 'new name'));

# test display_id

ok($vgf->display_id() eq $vgname);
