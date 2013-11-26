# Copyright [1999-2013] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


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
my $trv = $trva->fetch_by_dbID(9);

ok($trv->dbID() == 9);
ok(!defined($trv->cdna_start()));
ok(!defined($trv->cdna_end()));
ok(!defined($trv->translation_start()));
ok(!defined($trv->translation_end()));
ok(!defined($trv->pep_allele_string()));
ok($trv->consequence_type() eq 'DOWNSTREAM');
ok($trv->adaptor == $trva);


$trv = $trva->fetch_by_dbID(16);

ok($trv->dbID() == 16);
ok(!defined($trv->cdna_start()));
ok(!defined($trv->cdna_end()));
ok(!defined($trv->translation_start()));
ok(!defined($trv->translation_end()));
ok(!defined($trv->pep_allele_string()));
ok($trv->consequence_type() eq 'INTRONIC');
ok($trv->adaptor() == $trva);


# test fetch_all_by_VariationFeatures
my $slice = $sa->fetch_by_region('chromosome',20,30_600_000,31_000_000);
my $vf = $vfa->fetch_all_by_Slice($slice);
my @trvs = @{$trva->fetch_all_by_VariationFeatures($vf)};
ok(@trvs == 8);





