
use strict;
use warnings;

BEGIN { $| = 1;
	use Test;
	plan tests => 14;
}


use Bio::EnsEMBL::Test::TestUtils;


use Bio::EnsEMBL::Test::MultiTestDB;


our $verbose = 0;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();

my $vdb = $multi->get_DBAdaptor('variation');

my $aga = $vdb->get_AlleleGroupAdaptor();

ok($aga && $aga->isa('Bio::EnsEMBL::Variation::DBSQL::AlleleGroupAdaptor'));

# test fetch by dbID

my $ag = $aga->fetch_by_dbID(13146);

ok($ag->dbID() == 13146);
ok($ag->name() eq '376003');
my $vg = $ag->variation_group();
ok($vg->name() eq 'PERLEGEN:B003760');
ok($ag->source() eq 'dbSNP');
my @vars = @{$ag->get_all_Variations()};

ok($vars[0]->dbID() == 1292);
ok($vars[1]->dbID() == 1293);

my @alleles = @{$ag->get_all_alleles()};
ok($alleles[0] eq 'c');
ok($alleles[1] eq 'c');



# test fetch by name

$ag = $aga->fetch_by_name('ABDR-10');

ok($ag->dbID() == 10);
ok($ag->name() eq 'ABDR-10');


# test fetch all by Variation

my @ags = sort {$a->dbID <=> $b->dbID}
          @{$aga->fetch_all_by_VariationGroup($vg)};
ok(@ags == 5);
ok($ags[0]->dbID() == 13143);
ok(@{$ags[0]->get_all_Variations()} == 2);
