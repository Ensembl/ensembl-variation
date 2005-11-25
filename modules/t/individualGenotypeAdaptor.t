use strict;
use warnings;

BEGIN { $| = 1;
	use Test;
	plan tests => 8;
}


use Bio::EnsEMBL::Test::TestUtils;

use Bio::EnsEMBL::Test::MultiTestDB;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();

my $vdb = $multi->get_DBAdaptor('variation');

my $igty_adaptor = $vdb->get_IndividualGenotypeAdaptor();

ok($igty_adaptor->isa('Bio::EnsEMBL::Variation::DBSQL::CompressedGenotypeAdaptor'));

# # test fetch_all_by_individual
# my $ind_adaptor = $vdb->get_IndividualAdaptor();
# my $ind = shift @{$ind_adaptor->fetch_individual_by_synonym(1208)};

# my @igtys = sort {$a->variation->dbID() <=> $b->variation->dbID()}
#             @{$igty_adaptor->fetch_all_by_Individual($ind)};

# ok(@igtys == 17);
# ok($igtys[0]->variation()->name() eq 'rs193');
# ok($igtys[0]->allele1() eq 'C');
# ok($igtys[0]->allele2() eq 'T');
# ok($igtys[0]->individual()->dbID() == 1776);


# test fetch_all_by_Variation
my $variation_adaptor = $vdb->get_VariationAdaptor();
my $variation = $variation_adaptor->fetch_by_dbID(191);

my @igtys = ();

@igtys = sort {$a->individual->dbID() <=> $b->individual->dbID()}
            @{$igty_adaptor->fetch_all_by_Variation($variation)};

ok(@igtys == 50);
ok($igtys[0]->variation()->name() eq 'rs193');
ok($igtys[0]->allele1() eq 'C');
ok($igtys[0]->allele2() eq 'T');
ok($igtys[0]->individual()->name() eq 'NA17011');

#test get_all_Populations
my @pops = sort {$a->dbID() <=> $b->dbID()}
            @{$igtys[0]->individual->get_all_Populations()};
ok(@pops == 3);

ok($pops[0]->name eq 'TSC-CSHL:CEL_asian');
