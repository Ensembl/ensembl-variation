use strict;
use warnings;

BEGIN { $| = 1;
	use Test;
	plan tests => 13;
}


use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Variation::Variation;
use Bio::EnsEMBL::Variation::VariationFeature;




# test constructor

my $v = Bio::EnsEMBL::Variation::Variation->new(-name => 'rs2421',
                                                -source => 'dbSNP');


my $start = 100;
my $end   = 100;
my $strand = 1;
my $vname = $v->name();
my $map_weight = 1;
my $allele_str = 'A/T';

my $vf = Bio::EnsEMBL::Variation::VariationFeature->new
  (-start => 100,
   -end   => 100,
   -strand => 1,
   -variation_name => $vname,
   -map_weight => $map_weight,
   -allele_string => $allele_str,
   -variation => $v);


ok($vf->start()   == $start);
ok($vf->end()     == $end);
ok($vf->strand()  == $strand);
ok($vf->variation_name() eq $vname);
ok($vf->map_weight() == $map_weight);
ok($vf->allele_string() eq $allele_str);


# test getter/setters

my $v2 = Bio::EnsEMBL::Variation::Variation->new(-name => 'rs12311',
                                                 -source => 'dbSNP');

ok(test_getter_setter($vf, 'variation', $v2));
ok(test_getter_setter($vf, 'map_weight', 4));
ok(test_getter_setter($vf, 'allele_string', 'T|G'));
ok(test_getter_setter($vf, 'variation_name', 'rs21431'));


# test display_id

ok($vf->display_id eq $vname);

#test ambiguity code
ok($vf->ambig_code eq 'W');

#test variation class
ok($vf->var_class eq 'snp')
