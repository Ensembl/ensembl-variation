# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::Source;
our $verbose = 0;


my $omulti = Bio::EnsEMBL::Test::MultiTestDB->new('multi');
my $odb = $omulti->get_DBAdaptor('ontology');
Bio::EnsEMBL::Registry->add_db($omulti, 'ontology', $odb);

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();

my $vdb = $multi->get_DBAdaptor('variation');
my $db  = $multi->get_DBAdaptor('core');

$vdb->dnadb($db);

my $vf_ad  = $vdb->get_VariationFeatureAdaptor();
my $var_ad = $vdb->get_VariationAdaptor();
my $trv_ad = $vdb->get_TranscriptVariationAdaptor();
my $tr_ad  = $db->get_TranscriptAdaptor;
my $s_ad   = $db->get_SliceAdaptor();
my $g_ad = $db->get_GeneAdaptor();

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
ok( $trv->cdna_allele_string() eq 'C/G',                     "cdna_allele_string");
ok( $trv->affects_peptide() eq '1',                          "affects_peptide");
ok( $trv->hgvs_transcript()->{'G'} eq  'ENST00000470094.1:c.16C>G', "hgvs c");
ok( $trv->hgvs_protein()->{'G'}    eq  'ENSP00000434898.1:p.Ser6Cys', "hgvs p");

my $tvas = $trv->get_all_alternate_TranscriptVariationAlleles();
ok( $tvas->[0]->sift_prediction eq 'deleterious',            "sift prediction");


# test fetch_all_by_VariationFeatures
my $slice = $s_ad->fetch_by_region('chromosome',13,32953990,32954050);
my $vf = $vf_ad->fetch_all_by_Slice($slice);
my @trvs = @{$trv_ad->fetch_all_by_VariationFeatures($vf)};
ok(@trvs == 5 ,      "count by slice & variation features" );


is($trv_ad->fetch_by_dbID(938362444)->display_consequence, 'sequence_variant', 'empty consequence column');



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


my $tvs1 = $trv_ad->fetch_all_by_Transcripts_SO_terms([$transcript], ['sequence_variant']); 
ok(scalar @$tvs1 == 0, 'fetch_all_by_Transcripts_SO_terms');

my $tvs2 = $trv_ad->fetch_all_somatic_by_Transcripts_SO_terms([$transcript], ['sequence_variant']); 
ok(scalar @$tvs2 == 0, 'fetch_all_somatic_by_Transcripts_SO_terms');

my $tvs3 = $trv_ad->fetch_all_by_VariationFeatures_SO_terms([$new_vf], [$transcript], ['sequence_variant']); 
ok(scalar @$tvs3 == 0, 'fetch_all_by_VariationFeatures_SO_terms');

my $count4 = $trv_ad->count_all_by_VariationFeatures_SO_terms([$new_vf], [$transcript], ['sequence_variant']); 
ok($count4 == 0, 'count_all_by_VariationFeatures_SO_terms');

my $tvs5 = $trv_ad->fetch_all_somatic_by_Transcripts([$transcript]); 
ok(scalar @$tvs5 == 0, 'fetch_all_somatic_by_Transcripts');

my $translation_stable_id = $transcript->translation->stable_id;
ok($translation_stable_id eq 'ENSP00000434898', 'translation stable_id');

my $tvs6 = $trv_ad->fetch_all_by_translation_id($translation_stable_id);
ok(scalar @$tvs6 == 9, 'fetch_all_by_translation_id');

my $tvs7 = $trv_ad->fetch_all_somatic_by_translation_id($translation_stable_id);
ok(scalar @$tvs7 == 0, 'fetch_all_somatic_by_translation_id');

my $tvs8 = $trv_ad->fetch_all_by_translation_id_SO_terms($translation_stable_id, ['sequence_variant']);
ok(scalar @$tvs8 == 0, 'fetch_all_by_translation_id_SO_terms');

my $tvs9 = $trv_ad->fetch_all_somatic_by_translation_id_SO_terms($translation_stable_id, ['sequence_variant']);
ok(scalar @$tvs9 == 0, 'fetch_all_somatic_by_translation_id_SO_terms');




#### check HGVS shifting

my $source_name           = 'dbSNP';
my $source_version        = 138;
my $source_description    = 'Variants (including SNPs and indels) imported from dbSNP (mapped to GRCh38)';

my $source = Bio::EnsEMBL::Variation::Source->new
  (-name           => $source_name,
   -version        => $source_version,
   -description    => $source_description
);

## need a variation object 
my $v = Bio::EnsEMBL::Variation::Variation->new(-dbID => 12345,
                                                -name => 'rs2421',
                                                -source => $source);

my $chr   = 13;
my $start = 51519667;
my $end   = 51519668;
my $strand = 1;
my $vname = $v->name();
my $map_weight = 1;
my $allele_str = '-/G';
my $is_somatic = 0;


my $sl = $s_ad->fetch_by_region('chromosome',$chr);
$vf = Bio::EnsEMBL::Variation::VariationFeature->new
  (-seq_region_name => $chr,
   -start => $end,
   -end   => $start,
   -slice => $sl,
   -strand => $strand,
   -variation_name => $vname,
   -allele_string => $allele_str,
   -variation => $v,
   -source => $source,
   -is_somatic => $is_somatic,
   -adaptor     => $vf_ad
);

my $trans_name = 'ENST00000336617';
my $trans      = $tr_ad->fetch_by_stable_id( $trans_name);   

## test default - shifted
my $trans_vars = $vf->get_all_TranscriptVariations( [ $trans ] );

foreach my $trans_var (@{$trans_vars}){

  next unless $trans_var->transcript->stable_id() eq $trans_name;
  my $tvas_ts = $trans_var->get_all_alternate_TranscriptVariationAlleles();

  ok($tvas_ts->[0]->hgvs_transcript_reference() eq 'G', "hgvs_transcript_reference available without hgvs_transcript");
  ok($tvas_ts->[0]->hgvs_transcript() eq 'ENST00000336617.2:c.616+1dup', 'HGVS for shifted location');
  ok(scalar $tvas_ts->[0]->hgvs_offset() == 2, 'shifted offset');
}

## test non - shifted
$vf_ad->db->shift_hgvs_variants_3prime(0) ;
my $vf2 = Bio::EnsEMBL::Variation::VariationFeature->new
  (-seq_region_name => $chr,
   -start => $end,
   -end   => $start,
   -slice => $sl,
   -strand => $strand,
   -variation_name => $vname,
   -allele_string => $allele_str,
   -variation => $v,
   -source => $source,
   -is_somatic => $is_somatic,
   -adaptor     => $vf_ad
);

my $trans_vars_ns = $vf2->get_all_TranscriptVariations( [$trans] );

foreach my $trans_varns (@{$trans_vars_ns}){

  next unless $trans_varns->transcript->stable_id() eq $trans_name;
  my $tvas_ns = $trans_varns->get_all_alternate_TranscriptVariationAlleles(); 

  ok( $tvas_ns->[0]->hgvs_transcript() eq 'ENST00000336617.2:c.615_616insG', 'HGVS for non-shifted location' );
  ok(scalar $tvas_ns->[0]->hgvs_offset() == 0, 'non shifted offset');
}

#store
my $dbh = $vf_ad->dbc->db_handle;
my $sth = $dbh->prepare(qq/SELECT MAX(transcript_variation_id) FROM transcript_variation/); 
$sth->execute;
my ($max_tv_id_before) = $sth->fetchrow_array;
$sth->finish();

my $transcript_id_store = 'ENST00000470094';
my $transcript_store = $tr_ad->fetch_by_stable_id($transcript_id_store);
my $vf_store = $vf_ad->fetch_by_dbID(23700405);
my @vfs = ($vf_store);
foreach my $vf (@vfs) {
  my $tv = Bio::EnsEMBL::Variation::TranscriptVariation->new(
    -transcript     => $transcript_store,
    -variation_feature  => $vf,
    -adaptor      => $trv_ad,
    -disambiguate_single_nucleotide_alleles => 0,
    -no_transfer    => 1,
  );
  $trv_ad->store($tv);
}

sleep(10);

$sth->execute;
my ($max_tv_id_after) = $sth->fetchrow_array;
$sth->finish();

ok($max_tv_id_after == $max_tv_id_before + 1, 'get max transcript_variation_id');

my $tv_store = $trv_ad->fetch_by_dbID($max_tv_id_after);
ok($tv_store->display_consequence eq 'missense_variant', 'test store');
$dbh->do(qq{DELETE FROM variation_feature WHERE variation_feature_id=$max_tv_id_after;}) or die $dbh->errstr;

done_testing();

