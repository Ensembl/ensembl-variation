use strict;
use warnings;

BEGIN { $| = 1;
	use Test;
	plan tests => 8;
}


use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Variation::IndividualGenotype;
use Bio::EnsEMBL::Variation::Variation;
use Bio::EnsEMBL::Variation::Individual;


# test constructor
my $ind = Bio::EnsEMBL::Variation::Individual->new
  (-name => 'test individual',
   -description => 'This is a test individual',
   -gender => 'Male');

my $var = Bio::EnsEMBL::Variation::Variation->new
  (-name => 'rs123',
   -synonyms => {'dbSNP' => ['ss12', 'ss144']},
   -source => 'dbSNP');


my $a1 = 'A';
my $a2 = 'C';

my $ind_gtype = Bio::EnsEMBL::Variation::IndividualGenotype->new
  (-allele1 => $a1,
   -allele2 => $a2,
   -variation => $var,
   -individual => $ind);


ok($ind_gtype->individual()->name() eq $ind->name());
ok($ind_gtype->variation()->name()  eq $var->name());

ok($ind_gtype->allele1() eq $a1);
ok($ind_gtype->allele2() eq $a2);


# test getter/setters

ok(test_getter_setter($ind_gtype, 'allele1', 'TT'));
ok(test_getter_setter($ind_gtype, 'allele2', '-'));

my $ind2 = Bio::EnsEMBL::Variation::Individual->new
  (-name => 'test individual 2',
   -description => 'This is a second test individual',
   -gender => 'Female');

my $var2 = Bio::EnsEMBL::Variation::Variation->new
  (-name => 'rs332',
   -synonyms => {'dbSNP' => ['ss44', 'ss890']},
   -source => 'dbSNP');

ok(test_getter_setter($ind_gtype, 'individual', $ind2));
ok(test_getter_setter($ind_gtype, 'variation', $var2));
