use strict;
use warnings;

BEGIN { $| = 1;
	use Test;
	plan tests => 11;
}


use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Variation::PopulationGenotype;
use Bio::EnsEMBL::Variation::Variation;
use Bio::EnsEMBL::Variation::Population;


# test constructor
my $pop = Bio::EnsEMBL::Variation::Population->new
  (-name => 'test population',
   -description => 'This is a test individual',
   -size => 10000);

my $var = Bio::EnsEMBL::Variation::Variation->new
  (-name => 'rs123',
   -synonyms => {'dbSNP' => ['ss12', 'ss144']},
   -source => 'dbSNP');


my $a1 = 'A';
my $a2 = 'C';

my $dbID = 1;

my $freq = 0.76;

my $pop_gtype = Bio::EnsEMBL::Variation::PopulationGenotype->new
  (-dbID => $dbID,
   -allele1 => $a1,
   -allele2 => $a2,
   -variation => $var,
   -population => $pop,
   -frequency => $freq);


ok($pop_gtype->population()->name() eq $pop->name());
ok($pop_gtype->variation()->name()  eq $var->name());

ok($pop_gtype->allele1() eq $a1);
ok($pop_gtype->allele2() eq $a2);

ok($pop_gtype->dbID() == $dbID);

ok($pop_gtype->frequency() == $freq);

# test getter/setters

ok(test_getter_setter($pop_gtype, 'allele1', 'TT'));
ok(test_getter_setter($pop_gtype, 'allele2', '-'));

my $pop2 = Bio::EnsEMBL::Variation::Population->new
  (-name => 'test population 2',
   -description => 'This is a second test population',
   -size => 1000000);

my $var2 = Bio::EnsEMBL::Variation::Variation->new
  (-name => 'rs332',
   -synonyms => {'dbSNP' => ['ss44', 'ss890']},
   -source => 'dbSNP');

ok(test_getter_setter($pop_gtype, 'population', $pop2));
ok(test_getter_setter($pop_gtype, 'variation', $var2));

ok(test_getter_setter($pop_gtype, 'frequency', $freq));
