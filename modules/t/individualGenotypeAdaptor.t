use lib 't';
use strict;
use warnings;

BEGIN { $| = 1;
	use Test;
	plan tests => 12;
}


use TestUtils qw ( debug test_getter_setter count_rows);

use MultiTestDB;



my $multi = MultiTestDB->new();

my $vdb = $multi->get_DBAdaptor('variation');


my $igty_adaptor = $vdb->get_IndividualGenotypeAdaptor();

ok($igty_adaptor->isa('Bio::EnsEMBL::Variation::DBSQL::IndividualGenotypeAdaptor'));

# test fetch_by_dbID
my $igty = $igty_adaptor->fetch_by_dbID(5);

ok($igty->dbID() == 5);
ok($igty->allele1() eq 'A');
ok($igty->allele2() eq 'C');
ok($igty->variation()->name() eq 'rs241');
ok($igty->individual()->name() eq '1' && $igty->individual()->dbID() == 1208);


# test fetch_all_by_individual
my $ind_adaptor = $vdb->get_IndividualAdaptor();
my $ind = $ind_adaptor->fetch_by_dbID(1208);

my @igtys = sort {$a->dbID() <=> $b->dbID()}
            @{$igty_adaptor->fetch_all_by_Individual($ind)};

ok(@igtys == 16);
ok($igtys[0]->dbID() == 3);
ok($igtys[0]->variation()->name() eq 'rs221');
ok($igtys[0]->allele1() eq 'A');
ok($igtys[0]->allele2() eq 'G');
ok($igtys[0]->individual()->dbID() == 1208);


