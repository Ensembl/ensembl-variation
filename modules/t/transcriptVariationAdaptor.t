
use strict;
use warnings;

BEGIN { $| = 1;
	use Test;
	plan tests => 18;
}


use Bio::EnsEMBL::Test::TestUtils;


use Bio::EnsEMBL::Test::MultiTestDB;
our $verbose = 0;


my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();

my $vdb = $multi->get_DBAdaptor('variation');
my $db  = $multi->get_DBAdaptor('core');

$vdb->dnadb($db);

my $vfa = $vdb->get_VariationFeatureAdaptor();
my $va  = $vdb->get_VariationAdaptor();
my $trva = $vdb->get_TranscriptVariationAdaptor();
my $tra   = $db->get_TranscriptAdaptor;
my $sa = $db->get_SliceAdaptor();


ok($trva && $trva->isa('Bio::EnsEMBL::Variation::DBSQL::TranscriptVariationAdaptor'));



# test fetch_by_dbID
my $trv = $trva->fetch_by_dbID(1);

ok($trv->dbID() == 1);
ok(!defined($trv->cdna_start()));
ok(!defined($trv->cdna_end()));
ok(!defined($trv->translation_start()));
ok(!defined($trv->translation_end()));
ok(!defined($trv->pep_allele_string()));
ok($trv->consequence_type() eq 'UPSTREAM');
ok($trv->adaptor == $trva);


$trv = $trva->fetch_by_dbID(40);

ok($trv->dbID() == 40);
ok($trv->cdna_start() == 775);
ok($trv->cdna_end() == 775);
ok($trv->translation_start() == 255);
ok($trv->translation_end() == 255);
ok($trv->pep_allele_string() eq 'V/E');
ok($trv->consequence_type() eq 'NON_SYNONYMOUS_CODING');
ok($trv->adaptor() == $trva);


# test fetch_all_by_VariationFeatures
my $slice = $sa->fetch_by_region('chromosome',12,1,1e6);
my $vf = $vfa->fetch_all_by_Slice($slice);
my @trvs = @{$trva->fetch_all_by_VariationFeatures($vf)};
ok(@trvs == 2);





