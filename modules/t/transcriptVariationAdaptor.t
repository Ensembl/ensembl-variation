
use strict;
use warnings;

BEGIN { $| = 1;
	use Test;
	plan tests => 24;
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

ok($trva && $trva->isa('Bio::EnsEMBL::Variation::DBSQL::TranscriptVariationAdaptor'));



# test fetch_by_dbID

my $trv = $trva->fetch_by_dbID(1);

ok($trv->dbID() == 1);
# ok($trv->transcript->dbID() == 32371);
ok($trv->variation_feature->dbID() == 1039);
ok(!defined($trv->cdna_start()));
ok(!defined($trv->cdna_end()));
ok(!defined($trv->translation_start()));
ok(!defined($trv->translation_end()));
ok(!defined($trv->pep_allele_string()));
ok($trv->type() eq 'UPSTREAM');
ok($trv->adaptor == $trva);


$trv = $trva->fetch_by_dbID(40);

ok($trv->dbID() == 40);
# ok($trv->transcript->dbID() == 32829);
ok($trv->variation_feature->dbID() == 582);
ok($trv->cdna_start() == 775);
ok($trv->cdna_end() == 775);
ok($trv->translation_start() == 255);
ok($trv->translation_end() == 255);
ok($trv->pep_allele_string() eq 'V/E');
ok($trv->type() eq 'NON_SYNONYMOUS_CODING');
ok($trv->adaptor() == $trva);


# test fetch_all_by_VariationFeature
my $vf = $vfa->fetch_by_dbID(888);
my @trvs = @{$trva->fetch_all_by_VariationFeature($vf)};
ok(@trvs == 15);

# test fetch_all_by_Variation
my $v = $va->fetch_by_name('rs1002');
@trvs = @{$trva->fetch_all_by_Variation($v)};
ok(@trvs == 15);


# test fetch_all_by_Transcript

my $tr = $tra->fetch_by_stable_id('ENST00000252021');

@trvs = @{$trva->fetch_all_by_Transcript($tr)};
ok(@trvs == 1);
ok($trvs[0]->dbID() == 864);
ok($trvs[0]->transcript->dbID == $tr->dbID());




