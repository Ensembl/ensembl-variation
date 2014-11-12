# Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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
use Test::More;
use Data::Dumper;

use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Test::MultiTestDB;

our $verbose = 0;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');

my $vdb = $multi->get_DBAdaptor('variation');
my $db  = $multi->get_DBAdaptor('core');

$vdb->dnadb($db);

my $svfa = $vdb->get_StructuralVariationFeatureAdaptor();
my $sva  = $vdb->get_StructuralVariationAdaptor();

my $dbID = 4509635;
my $outer_start = 7803891;
my $inner_start = 7805991;
my $inner_end = 7823440;
my $outer_end = 7825340;
my $var_name = 'esv93078';
my $source_name = 'DGVa';
my $study_name  = 'estd59';
my $chr = '8';

ok($svfa && $svfa->isa('Bio::EnsEMBL::Variation::DBSQL::StructuralVariationFeatureAdaptor'));

my $svf = $svfa->fetch_by_dbID($dbID);
my $source = $svf->source_object;
my $slice  = $svf->slice;

ok($svf->outer_start() == $outer_start, "svf_id -> outer_start");
ok($svf->start()       == $outer_start, "svf_id -> start");
ok($svf->inner_start() == $inner_start, "svf_id -> inner_start");
ok($svf->inner_end()   == $inner_end,   "svf_id -> inner_end");
ok($svf->end()         == $outer_end,   "svf_id -> end") ;
ok($svf->outer_end()   == $outer_end,   "svf_id -> outer_end");
ok($svf->strand()      == 1,            "svf_id -> strand");
ok($svf->variation_name() eq $var_name, "svf_id -> varname" );
ok($source->name() eq $source_name,     "svf_id -> source" );
ok($svf->study->name() eq $study_name,  "svf_id -> study" );
ok($slice->seq_region_name() eq $chr,   "svf_id -> slice name");



# test fetch_all_by_Variation

my $sv = $sva->fetch_by_dbID(3506221);
my $svfs = $svfa->fetch_all_by_StructuralVariation($sv);
ok(@$svfs == 1,                         "sv -> vf count ");
$svf = $svfs->[0];
$source = $svf->source_object;
$slice  = $svf->slice;

ok($svf->dbID() == $dbID,               "sv -> svf id");
ok($slice->seq_region_name() eq $chr,   "sv -> slice name ");
ok($svf->outer_start() == $outer_start, "sv -> outer_start");
ok($svf->start()       == $outer_start, "sv -> start");
ok($svf->inner_start() == $inner_start, "sv -> inner_start");
ok($svf->inner_end()   == $inner_end,   "sv -> inner_end");
ok($svf->end()         == $outer_end,   "sv -> end") ;
ok($svf->outer_end()   == $outer_end,   "sv -> outer_end");
ok($svf->strand()      == 1,            "sv -> strand");
ok($svf->variation_name() eq $var_name, "sv -> name");
ok($source->name() eq $source_name,     "sv -> source" );
ok($svf->study->name() eq $study_name,  "sv -> study" );

done_testing();

