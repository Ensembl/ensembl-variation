# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2020] EMBL-European Bioinformatics Institute
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
use Test::Exception;
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
my $ofdb  = $multi->get_DBAdaptor('otherfeatures');

$vdb->dnadb($db);

my $vf_ad  = $vdb->get_VariationFeatureAdaptor();
my $var_ad = $vdb->get_VariationAdaptor();
my $trv_ad = $vdb->get_TranscriptVariationAdaptor();
my $tr_ad  = $db->get_TranscriptAdaptor;
my $s_ad   = $db->get_SliceAdaptor();
my $g_ad = $db->get_GeneAdaptor();
my $tr_ad_of = $ofdb->get_TranscriptAdaptor();

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
ok( $trv->cdna_start_unshifted() eq 16,                     "cdna_start_unshifted");
ok( $trv->cdna_end_unshifted()   eq 16,                     "cdna_end_unshifted");
ok( $trv->translation_start_unshifted() eq 6,               "translation_start_unshifted");
ok( $trv->translation_end_unshifted() eq 6,                 "translation_end_unshifted");
ok( $trv->pep_allele_string() eq 'S/C',                     "pep_allele_string");
ok( $trv->consequence_type()->[0] eq 'missense_variant',    "consequence");
ok( $trv->variation_feature->variation_name() eq "rs80359157", "variation name ");
ok( $trv->transcript->stable_id() eq "ENST00000470094",      "feature id" );
ok( $trv->cdna_allele_string() eq 'C/G',                     "cdna_allele_string");
ok( $trv->affects_peptide() eq '1',                          "affects_peptide");
ok( $trv->hgvs_transcript()->{'G'} eq  'ENST00000470094.1:c.16C>G', "hgvs c");
ok( $trv->hgvs_protein()->{'G'}    eq  'ENSP00000434898.1:p.Ser6Cys', "hgvs p");
my $tvas = $trv->get_all_alternate_TranscriptVariationAlleles();
my $tva = $tvas->[0];
ok( $tva->sift_prediction eq 'deleterious',            "sift prediction");
ok($tva->dbnsfp_meta_lr_prediction eq 'tolerated', "dbnsfp meta_lr prediction");
ok($tva->dbnsfp_meta_lr_score == 0.459, "dbnsfp meta_lr score");
ok($tva->dbnsfp_mutation_assessor_prediction eq 'medium', "dbnsfp mutation_assessor prediction");
ok($tva->dbnsfp_mutation_assessor_score == 0.767, "dbnsfp mutation_assessor score");
ok($tva->dbnsfp_revel_prediction eq 'likely benign', "dbnsfp revel prediction");
ok($tva->dbnsfp_revel_score == 0.315, "dbnsfp revel score");
ok($tva->cadd_prediction eq 'likely benign', "cadd prediction");
ok($tva->cadd_score == 23, "cadd score");

$tva->clear_shifting_variables();
ok(!defined($trv->{cds_start}), 'cleared cds_start');
ok(!defined($trv->{cds_end}), 'cleared cds_end');
$tva->{shift_hash}->{shift_length} = 0;
$tva->clear_shifting_variables();

ok(defined($trv->{cds_start}), 'redefined cds_start');
ok(defined($trv->{cds_end}), 'redefined cds_end');
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

my $count10 = $trv_ad->count_all_by_Transcript( $transcript );
ok( $count10 == 10, 'count_all_by_Transcript');

throws_ok { $trv_ad->count_all_by_Transcript(); } qr/Bio::EnsEMBL::Transcript arg expected/, "count missing transcript";

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
  ok($tvas_ts->[0]->hgvs_intron_start_offset == 1, 'hgvs_intron_start_offset');
  ok($tvas_ts->[0]->hgvs_exon_start_coordinate == 616, 'hgvs_exon_start_coordinate');
  ok($tvas_ts->[0]->hgvs_intron_end_offset == 1, 'hgvs_intron_end_offset');
  ok($tvas_ts->[0]->hgvs_exon_end_coordinate == 616, 'hgvs_exon_end_coordinate');
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

#Test shifting without db adaptor

$vf_ad->db->shift_hgvs_variants_3prime(0) ;
my $vf3 = Bio::EnsEMBL::Variation::VariationFeature->new
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

my $trans_vars_ns2 = $vf3->get_all_TranscriptVariations( [$trans] );

delete($vf3->{adaptor});
$Bio::EnsEMBL::Variation::DBSQL::DBAdaptor::DEFAULT_SHIFT_HGVS_VARIANTS_3PRIME = 0;

foreach my $trans_varns (@{$trans_vars_ns2}){

  next unless $trans_varns->transcript->stable_id() eq $trans_name;
  my $tvas_ns = $trans_varns->get_all_alternate_TranscriptVariationAlleles();

  ok( $tvas_ns->[0]->hgvs_transcript() eq 'ENST00000336617.2:c.615_616insG', 'HGVS for non-shifted location without adaptor' );
  ok(scalar $tvas_ns->[0]->hgvs_offset() == 0, 'non shifted offset without adaptor');
}


#store

#make a backup of the current transcript_variation table before the store test
$multi->save('variation', 'transcript_variation');

my $dbh = $vf_ad->dbc->db_handle;
my $sth = $dbh->prepare(qq/SELECT MAX(transcript_variation_id) FROM transcript_variation/);
$sth->execute;
my ($max_tv_id_before) = $sth->fetchrow_array;
$sth->finish();

my $transcript_id_store = 'ENST00000470094';
my $transcript_store = $tr_ad->fetch_by_stable_id($transcript_id_store);
my $vf_store = $vf_ad->fetch_by_dbID(23700405);
my $tv = Bio::EnsEMBL::Variation::TranscriptVariation->new(
    -transcript     => $transcript_store,
    -variation_feature  => $vf_store,
    -adaptor      => $trv_ad,
    -disambiguate_single_nucleotide_alleles => 0,
    -no_transfer    => 1,
);
$trv_ad->store($tv);

sleep(10);

$sth->execute;
my ($max_tv_id_after) = $sth->fetchrow_array;
$sth->finish();

# When records are deleted, the auto_increment of the table is not changed.
# When running with database intact, the max_tv_id_after should be
# greater than before but not necessarily one greater
ok($max_tv_id_after > $max_tv_id_before, 'get max transcript_variation_id');

my $tv_store = $trv_ad->fetch_by_dbID($max_tv_id_after);
ok($tv_store->display_consequence eq 'missense_variant', 'test store');

# restore the transcript_variation table from before store test
$multi->restore('variation', 'transcript_variation');
# test most_severe_OverlapConsequence for consequences with the same rank
$stable_id = 'ENST00000470094';
$transcript = $tr_ad->fetch_by_stable_id($stable_id);

my $var_msc = $var_ad->fetch_by_name('rs200814465');
my $vfs_msc = $vf_ad->fetch_all_by_Variation($var_msc);
my $new_vf_msc = $vfs_msc->[0];
my $trvar_msc = Bio::EnsEMBL::Variation::TranscriptVariation->new
  (-variation_feature => $new_vf_msc,
   -transcript        => $transcript,
);

# The expected consequences are:
#   rank: 12  SO_term: missense_variant
#   rank: 22  SO_term: NMD_transcript_variant
# For the test change to the same rank, different description
for my $allele (@{$trvar_msc->get_all_alternate_BaseVariationFeatureOverlapAlleles}) {
    for my $cons (@{$allele->get_all_OverlapConsequences }) {
      if ($cons->SO_term eq 'missense_variant') {
        $cons->SO_term('test_consequence_z');
        $cons->rank(3);
      } elsif ($cons->SO_term eq 'NMD_transcript_variant') {
        $cons->SO_term ('test_consequence_a');
        $cons->rank(3);
      }
    }
}
my $msc_2_expected = 'test_consequence_a';
my $msc_2 = $trvar_msc->most_severe_OverlapConsequence();
is($msc_2->SO_term, $msc_2_expected, 'tv - most_severe_OverlapConsequence - same rank');

## RefSeq Mismatch Testing
my $tr = $tr_ad_of->fetch_by_stable_id('NM_001270408.1');
my $sl_refseq = $s_ad->fetch_by_region('chromosome', 21);
$vf_ad->db->shift_hgvs_variants_3prime(0);

my $vf_refseq = Bio::EnsEMBL::Variation::VariationFeature->new
  (-seq_region_name => 21,
   -start => 25700000,
   -end   => 25700000,
   -slice => $sl_refseq,
   -strand => 1,
   -variation_name => 'refseq_test',
   -allele_string => 'G/A',
   -adaptor     => $vf_ad
);

$tv = Bio::EnsEMBL::Variation::TranscriptVariation->new(
    -transcript     => $tr,
    -variation_feature  => $vf_refseq,
    -adaptor      => $trv_ad,
    -no_shift     => 1,
);

## Setting CDS coordinates allows us to set up niche case where mismatches are calculated
my @tvas = @{ $tv->get_all_alternate_TranscriptVariationAlleles };

## The following predicate is normally set in _bvfo_preds, however due to the limited quantity of data within the test database,
## calling the method is insufficient right now. Rather than rewriting chunks of the test db, I've added in the predicate value here
$tvas[0]->{pre_consequence_predicates}->{exon} = 1;

$tv->cds_start(1234);
$tv->cds_end(1233);
$tv->cdna_start(1000);
$tv->transcript->{cdna_coding_start} = 234;

## Coordinate within HGVS matches cds_start as these values take mismatch into account 
ok($tv->hgvs_transcript->{A} eq 'NM_001270408.1:c.1234N>A', 'Refseq HGVS mismatch calculated');

## Misalignment offset calculated from transcript edits recognises insertion of 4BP
my @attribs = @{$tr->get_all_Attributes()};
my @edit_attrs = grep {$_->code =~ /^_rna_edit/} @attribs;
ok($tvas[0]->get_misalignment_offset(\@edit_attrs) == 4, 'Misalignment offset');


# Testing for transcript_variation.hgvs_genomic 'NULL' (string)
# and NULL (value)
#
# An existing database record is used for testing
# A backup of transcript_variation done at start of test
# and restored at end of test.
{
  # Backup the table before test
  $multi->save('variation', 'transcript_variation');

  # Data used for test
  # transcript_variation_id: 92216616
  #          variation_name: rs7949754
  #           allele_string: G/A
  #    variation_feature_id: 4806434
  # transcript_variation_id: 92216616
  #            hgvs_genomic: 11:g.66328199G>A
  my $dbID = 92216616;
  my $variant = 'rs7949754';
  my $expected_hgvs_genomic = 'NC_000011.9:g.66328199G>A';
  my $dbh = $trv_ad->dbc->db_handle();

  # Test setting transcript_variation.hgvs_genomic
  # to 'NULL' (the string)
  my $sql = qq{
    UPDATE transcript_variation
    SET hgvs_genomic = 'NULL'
    WHERE transcript_variation_id = $dbID
  };

  my $count = $dbh->do($sql);
  ok($count == 1 , "Update to 'NULL' string for testing");

  # Fetch by dbID
  my $tv_by_dbID = $trv_ad->fetch_by_dbID($dbID);
  my $hgvs_db = $tv_by_dbID->{alt_alleles}->[0]->{hgvs_genomic};
  is($hgvs_db, undef, "hgvs_genomic set to undefined when 'NULL'");

  # Fetch all TranscriptVariations related to a VariationFeature
  my $var = $var_ad->fetch_by_name($variant);
  my $vf = $vf_ad->fetch_all_by_Variation($var_ad->fetch_by_name($variant))->[0];
  my $tv = @{$trv_ad->fetch_all_by_VariationFeatures([$vf])}[0];
  my $hgvs_genomic = $tv->hgvs_genomic->{'A'};
  is($hgvs_genomic, $expected_hgvs_genomic, "Lookup hgvs_genomic when 'NULL' in TV");

  # Testing setting transcript_variation.hgvs_genomic to NULL (the value)
  $sql = qq{
    UPDATE transcript_variation
    SET hgvs_genomic = NULL
    WHERE transcript_variation_id = $dbID
  };

  $count = $dbh->do($sql);
  ok($count == 1 , "Update to NULL for testing");

  $tv_by_dbID = $trv_ad->fetch_by_dbID($dbID);
  $hgvs_db = $tv_by_dbID->{alt_alleles}->[0]->{hgvs_genomic};
  is($hgvs_db, undef, 'hgvs_genomic set to undefined when NULL');

  # Fetch all TranscriptVariations related to a VariationFeature
  my $tv_2 = @{$trv_ad->fetch_all_by_VariationFeatures([$vf])}[0];
  $hgvs_genomic = $tv_2->hgvs_genomic->{'A'};
  is($hgvs_genomic, $expected_hgvs_genomic, "Lookup hgvs_genomic when NULL in TV");

  # Restore table to what it was at start of test
  $multi->restore('variation', 'transcript_variation');
}

done_testing();
