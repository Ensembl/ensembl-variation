use lib 't';

use strict;
use warnings;

BEGIN { $| = 1;
	use Test;
	plan tests => 10;
}


use TestUtils qw ( debug test_getter_setter count_rows);


use MultiTestDB;


our $verbose = 0;

my $multi = MultiTestDB->new();

my $vdb = $multi->get_DBAdaptor('variation');


my $pa = $vdb->get_PopulationAdaptor();
ok($pa && $pa->isa('Bio::EnsEMBL::Variation::DBSQL::PopulationAdaptor'));

# test fetch by dbID

my $pop = $pa->fetch_by_dbID(7);

ok($pop->name() eq 'PACIFIC');
ok($pop->dbID() == 7);
debug("Descrip = '". $pop->description."'");
ok($pop->description() eq ' Samples from Australia, New Zealand, Central and Southern Pacific Islands, Southeast Asian Peninsular/Island Nations.');


# test fetch by name

$pop = $pa->fetch_by_name('NUSPAE:Chinese_HDL');
ok($pop->name() eq 'NUSPAE:Chinese_HDL');
ok($pop->dbID() == 97);

# test fetch_by_super_Population
$pop = $pa->fetch_by_name('PACIFIC');
my $pops = $pa->fetch_all_by_super_Population($pop);
ok(@$pops == 25);
ok($pops->[0]->dbID());

# test fetch_by_sub_Population
$pop = $pa->fetch_by_name('WRAYLAB:PNG');
$pops = $pa->fetch_all_by_sub_Population($pop);
ok(@$pops == 1);
ok($pops->[0]->name() eq 'PACIFIC');
