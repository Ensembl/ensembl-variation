use strict;
use warnings;

BEGIN { $| = 1;
	use Test;
	plan tests => 20;
}


use Bio::EnsEMBL::Test::TestUtils;

use Bio::EnsEMBL::Test::MultiTestDB;



my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();

my $vdb = $multi->get_DBAdaptor('variation');


my $pgty_adaptor = $vdb->get_PopulationGenotypeAdaptor();

ok($pgty_adaptor->isa('Bio::EnsEMBL::Variation::DBSQL::PopulationGenotypeAdaptor'));

# test fetch_by_dbID
my $pgty = $pgty_adaptor->fetch_by_dbID(5);

ok($pgty->dbID() == 5);
ok($pgty->allele1() eq 'C');
ok($pgty->allele2() eq 'G');
ok($pgty->variation()->name() eq 'rs2872');
ok($pgty->population()->name() eq 'AFFY:Caucasian');
ok($pgty->frequency() == 0.454545);


# test fetch_all_by_Population
my $pop_adaptor = $vdb->get_PopulationAdaptor();
my $pop = $pop_adaptor->fetch_by_dbID(429);

my @igtys = sort {$a->dbID() <=> $b->dbID()}
            @{$pgty_adaptor->fetch_all_by_Population($pop)};

ok(@igtys == 3);
ok($igtys[0]->dbID() == 1);
ok($igtys[0]->variation()->name() eq 'rs2872');
ok($igtys[0]->allele1() eq 'C');
ok($igtys[0]->allele2() eq 'C');
ok($igtys[0]->frequency() == 0.666667);
ok($igtys[0]->population()->dbID() == 429);


#test fetch_all_by_Variation

my $variation_adaptor = $vdb->get_VariationAdaptor();
my $variation = $variation_adaptor->fetch_by_dbID(2863);

@igtys = ();

@igtys = sort {$a->dbID() <=> $b->dbID()}
            @{$pgty_adaptor->fetch_all_by_Variation($variation)};

ok(@igtys == 12);
ok($igtys[0]->dbID() == 1);
ok($igtys[0]->population()->name() eq 'AFFY:AfAm');
ok($igtys[0]->allele1() eq 'C');
ok($igtys[0]->allele2() eq 'C');
ok($igtys[0]->frequency() == 0.666667);

