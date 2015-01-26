# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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


my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();

my $vdb = $multi->get_DBAdaptor('variation');
my $db  = $multi->get_DBAdaptor('core');

$vdb->dnadb($db);

my $vf_ad  = $vdb->get_VariationFeatureAdaptor();
my $var_ad = $vdb->get_VariationAdaptor();
my $trv_ad = $vdb->get_TranscriptVariationAdaptor();
my $tr_ad  = $db->get_TranscriptAdaptor;
my $s_ad   = $db->get_SliceAdaptor();


ok($trv_ad && $trv_ad->isa('Bio::EnsEMBL::Variation::DBSQL::TranscriptVariationAdaptor'));

## trans: ENST00000470094


#use API tutorial example:
my $stable_id = 'ENST00000470094'; #this is the stable_id of a human transcript

my $transcript = $tr_ad->fetch_by_stable_id($stable_id); #get the Transcript object

my $trvs = $trv_ad->fetch_all_by_Transcripts([$transcript]); #get ALL effects of Variations in the Transcript


ok(@$trvs == 9 , "transcriptvariation count") ;


my $trv = $trv_ad->fetch_by_dbID(49961554);

ok( $trv->dbID() == 49961554,                               "dbID");
ok( $trv->cdna_start() eq 16,                               "cdna_start");
ok( $trv->cdna_end()   eq 16,                               "cdna_end");
ok( $trv->translation_start() eq 6,                         "translation_start");
ok( $trv->translation_end() eq 6,                           "translation_end");
ok( $trv->pep_allele_string() eq 'S/C',                     "pep_allele_string");
ok( $trv->consequence_type()->[0] eq 'missense_variant',    "consequence");
ok( $trv->variation_feature->variation_name() eq "rs80359157", "variation name ");
ok( $trv->transcript->stable_id() eq "ENST00000470094",      "feature id" );

my $tvas = $trv->get_all_alternate_TranscriptVariationAlleles();
ok( $tvas->[0]->sift_prediction eq 'deleterious',            "sift prediction");




# test fetch_all_by_VariationFeatures
my $slice = $s_ad->fetch_by_region('chromosome',13,32953990,32954050);
my $vf = $vf_ad->fetch_all_by_Slice($slice);
my @trvs = @{$trv_ad->fetch_all_by_VariationFeatures($vf)};
ok(@trvs == 5 ,      "count by slice & variation features" );




# test constructor

my $var = $var_ad->fetch_by_name("rs200814465");

my $vfs = $vf_ad->fetch_all_by_Variation($var);
my $new_vf = $vfs->[0];
my $pep_allele = 'Q/H';

my $cdna_start = 68;
my $cdna_end   = 68;
my $tl_start   = 23;
my $tl_end     = 23;
my $consequence_type = 'missense_variant';


my $trvar = Bio::EnsEMBL::Variation::TranscriptVariation->new
  (-variation_feature => $new_vf,
   -transcript        => $transcript,
);



ok($trvar->variation_feature() == $new_vf,      "constructor - VF");
ok($trvar->pep_allele_string() eq $pep_allele,  "constructor - pep_allele_string");
ok($trvar->cdna_start() == $cdna_start,         "constructor - cdna_start");
ok($trvar->cdna_end() == $cdna_end,             "constructor - cdna_end");
ok($trvar->translation_start() == $tl_start,    "constructor - translation start");
ok($trvar->translation_end() == $tl_end,        "constructor - translation end");
ok($trvar->consequence_type()->[0] eq $consequence_type, "constructor - consequence_type");



# test getter/setters
my $tr_new = Bio::EnsEMBL::Transcript->new();
ok(test_getter_setter($trvar, 'transcript', $tr_new));




ok(test_getter_setter($trvar, 'pep_allele_string', $pep_allele));

ok(test_getter_setter($trvar, 'cdna_start', 1));
ok(test_getter_setter($trvar, 'cdna_end', 12));

ok(test_getter_setter($trvar, 'translation_start', 4));
ok(test_getter_setter($trvar, 'translation_end', 10));


done_testing();

