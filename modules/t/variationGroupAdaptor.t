
use strict;
use warnings;

BEGIN { $| = 1;
	use Test;
	plan tests => 15;
}


use Bio::EnsEMBL::Test::TestUtils;


use Bio::EnsEMBL::Test::MultiTestDB;


our $verbose = 0;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();

my $vdb = $multi->get_DBAdaptor('variation');

my $vga = $vdb->get_VariationGroupAdaptor();

ok($vga && $vga->isa('Bio::EnsEMBL::Variation::DBSQL::VariationGroupAdaptor'));

# test fetch by dbID

my $vg = $vga->fetch_by_dbID(3761);

ok($vg->dbID == 3761);
ok($vg->adaptor == $vga);
ok($vg->name eq 'PERLEGEN:B003760');
ok($vg->source eq 'dbSNP');
ok($vg->type eq 'haplotype');

my @vars = sort {$a->dbID() <=> $b->dbID} @{$vg->get_all_Variations()};

ok(@vars == 2);
ok($vars[0]->dbID() == 1292);
ok($vars[1]->dbID() == 1293);


# test fetch by name

$vg = $vga->fetch_by_name('DBMHC:ABDR');
ok($vg->dbID() == 1);
ok($vg->name() eq 'DBMHC:ABDR');
ok(!@{$vg->get_all_Variations});

# test fetch_all_by_Variation
my $vgs = $vga->fetch_all_by_Variation($vars[0]);

ok(@$vgs == 1);
ok($vgs->[0]->dbID() == 3761);
ok($vgs->[0]->name() eq 'PERLEGEN:B003760');
