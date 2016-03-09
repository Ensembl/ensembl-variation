# Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

use Bio::EnsEMBL::Variation::TranscriptHaplotypeContainer;

use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Test::MultiTestDB;


our $verbose = 0;


## SETUP
########

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');

my $vdb = $multi->get_DBAdaptor('variation');
my $cdb = $multi->get_DBAdaptor('core');

# set the VCFCollection config
my $dir = $multi->curr_dir();
ok($vdb->vcf_config_file($dir.'/th_vcf_config.json') eq $dir.'/th_vcf_config.json', "DBAdaptor vcf_config_file");
my $vca = $vdb->get_VCFCollectionAdaptor();
my $vcf_coll = $vca->fetch_all->[0];
my $temp = $vcf_coll->filename_template();
$temp =~ s/###t\-root###/$dir/;
$vcf_coll->filename_template($temp);

# get transcript
my $ta = $cdb->get_TranscriptAdaptor;
my $tr = $ta->fetch_by_stable_id('ENST00000502692');

my $tha = $vdb->get_TranscriptHaplotypeAdaptor();
my $c = $tha->get_TranscriptHaplotypeContainer_by_Transcript($tr);

my ($s) = grep {$_->name =~ /NA18499/} @{$vcf_coll->get_all_Samples};
ok($s && $s->isa('Bio::EnsEMBL::Variation::Sample'), "get_all_Samples");


## TESTS
########

my ($h, @h, $counts, $exp, $diffs);

# create container
ok($c->isa('Bio::EnsEMBL::Variation::TranscriptHaplotypeContainer'), "create container");

ok($c->db->isa('Bio::EnsEMBL::Variation::DBSQL::DBAdaptor'), "db method");



# get methods
$h = $c->_get_TranscriptHaplotype_by_hex('25f6ab7d0476132b731c630d743d2868');
ok($h->isa('Bio::EnsEMBL::Variation::CDSHaplotype'), "get by hex - isa");
is($h->name, 'ENST00000502692:1697G>A,1818A>G,1858T>C,2011T>C,2033A>C,2456A>G', "get by hex - name");

$h = $c->get_TranscriptHaplotype_by_name('ENST00000502692:1697G>A,1818A>G,1858T>C,2011T>C,2033A>C,2456A>G');
ok($h->isa('Bio::EnsEMBL::Variation::CDSHaplotype'), "get by name - isa");
is($h->_hex, '25f6ab7d0476132b731c630d743d2868', "get by name - hex");

@h = @{$c->get_all_TranscriptHaplotypes};
is(scalar @h, 116, "get all TranscriptHaplotypes - count");

@h = @{$c->get_all_CDSHaplotypes};
is(scalar @h, 75, "get all CDSHaplotypes - count");

@h = @{$c->get_all_ProteinHaplotypes};
is(scalar @h, 41, "get all ProteinHaplotypes - count");

@h = @{$c->get_all_TranscriptHaplotypes_by_Sample($s)};
is(scalar @h, 4, "get_all_TranscriptHaplotypes_by_Sample - count");

@h = @{$c->get_all_CDSHaplotypes_by_Sample($s)};
is(scalar @h, 2, "get_all_CDSHaplotypes_by_Sample - count");

@h = @{$c->get_all_ProteinHaplotypes_by_Sample($s)};
is(scalar @h, 2, "get_all_ProteinHaplotypes_by_Sample - count");

@h = @{$c->get_all_most_frequent_CDSHaplotypes};
ok(scalar @h == 1, "get_all_most_frequent_CDSHaplotypes");
is($h[0]->name, 'ENST00000502692:1697G>A,1858T>C,2011T>C,2456A>G', "get_all_most_frequent_CDSHaplotypes - name");

@h = @{$c->get_all_most_frequent_ProteinHaplotypes};
ok(scalar @h == 1, "get_all_most_frequent_ProteinHaplotypes");
is($h[0]->name, 'ENSP00000422007:566R>Q,620*>R,671C>R,819Q>R', "get_all_most_frequent_ProteinHaplotypes - name");

@h = @{$c->get_all_Populations};
is(scalar @h, 31, "get_all_Populations - count");



# counts etc
is($c->total_haplotype_count, 5008, "total haplotype count");

$counts = $c->total_population_counts();
$exp = {
  '1000GENOMES:phase_3_YRI' => 216,
  '1000GENOMES:phase_3_GIH' => 206,
  '1000GENOMES:phase_3_GWD' => 226,
  '1000GENOMES:phase_3_PUR' => 208,
  '1000GENOMES:phase_3_MSL' => 170,
  '1000GENOMES:phase_3_ITU' => 204,
  '1000GENOMES:phase_3_ESN' => 198,
  '1000GENOMES:phase_3_ACB' => 192,
  '1000GENOMES:phase_3_PEL' => 170,
  '1000GENOMES:phase_3_CEU' => 198,
  '1000GENOMES:phase_3_EAS' => 1008,
  '1000GENOMES:phase_3_FIN' => 198,
  '1000GENOMES:phase_3_CDX' => 186,
  '1000GENOMES:phase_3_AFR' => 1322,
  '1000GENOMES:phase_3_PJL' => 192,
  '1000GENOMES:phase_3_STU' => 204,
  '1000GENOMES:phase_3_LWK' => 198,
  '1000GENOMES:phase_3_IBS' => 214,
  '1000GENOMES:phase_3_CLM' => 188,
  '1000GENOMES:phase_3_EUR' => 1006,
  '1000GENOMES:phase_3_SAS' => 978,
  '1000GENOMES:phase_3_ASW' => 122,
  '1000GENOMES:phase_3_CHS' => 210,
  '1000GENOMES:phase_3_TSI' => 214,
  '1000GENOMES:phase_3_MXL' => 128,
  '1000GENOMES:phase_3_GBR' => 182,
  '1000GENOMES:phase_3_CHB' => 206,
  '1000GENOMES:phase_3_AMR' => 694,
  '1000GENOMES:phase_3_BEB' => 172,
  '1000GENOMES:phase_3_JPT' => 208,
  '1000GENOMES:phase_3_KHV' => 198
};
is_deeply($counts, $exp, "total population counts");


## TranscriptHaplotype tests

@h = @{$c->get_all_ProteinHaplotypes_by_Sample($s)};
$h = $h[1];

is($h->container, $c, 'container');
is($h->type, 'protein', 'type');
is($h->seq, 'MGIGLCVPPARCWKRLTRSSRGQPLERPAATIIIGLSIVERSKQTQRGPVDFSRSHSAPAAEQECKPRISAARSTPAPPSPAVPGFRTTPCQTFTAWCNSHLRKAGTQIENIEEDFRNGLKLMLLLEVISGERLPRPDKGKMRFHKIANVNKALDFIASKGVKLVSIGAEEIVDGNLKMTLGMIWTIILRFAIQDISVEETSAKEGLLLWCQRKTAPYRNVNVQNFHTSWKDGLALCALIHRHRPDLIDYAKLRKDDPIGNLNTAFEVAEKYLDIPKMLDAEDIVNTPKPDEKAIMTYVSCFYHAFAGAEQAETAANRICKVLAVNQENEKLMEEYEKLASELLEWIRRTVPWLENRVGEPSMSAMQRKLEDFRDYRRLHKPPRIQEKCQLEINFNTLQTKLRLSHRPAFMPSEGKLVSDIANAWRGLEQVEKGYEDWLLSEIRRLQRLQHLAEKFRQKASLHEAWTRGKEEMLSQRDYDSALLQEVRALLRRHEAFESDLAAHQDRVEHIAALAQELNELDYHEAASVNSRCQAICDQWDNLGTLTQKRRDALERMEKLLETIDQLQLEFARRAAPFNNWLDGAVEDLQDVWLVHSVEETQSLLTAHDQFKATLPEADRERGAIMGIQGEIQKICQTYGLRPCSTNPYITLSPQDINTKWDMVRKLVPSRDQTLQEALARQQVNERLRRQFAAQANAIGPWIQAKVEEVGRLAAGLAGSLEEQMAGLRQQEQNIINYKTNIDRLEGDHQLLQESLVFDNKHTVYSMEHIRVGWEQLLTSIARTINEVENQVLTRDAKGLSQEQLNEFRASFNHFDRKRNGMMEPDDFRACLISMGYDLGEVEFARIMTMVDPNAAGVVTFQAFIDFMTRETAETDTTEQVVASFKILAGDKNYITPELLRRELPAKQAEYCIRRMVPYKGSGAPAGALDYVAFSSALYGESDL*', 'seq');
is_deeply($h->get_all_flags, ['deleterious_sift_or_polyphen','stop_change'], 'get_all_flags');
is($h->has_indel, 0, 'has_indel');
is($h->has_deleterious_sift_or_polyphen, 1, 'has_deleterious_sift_or_polyphen');
is($h->has_stop_change, 1, 'has_stop_change');
is($h->_hex, 'dff67d1fb7b3812a10cf152d2d9528f7', '_hex');
is_deeply($h->_other_hexes, ['fbb87de3331ea4225df6847e8335d622'], '_other_hexes');
is($h->transcript, $c->transcript, 'transcript');
is($h->count, 47, 'count');
is($h->frequency, 0.00938498402555911, 'frequency');

is_deeply($h->get_all_population_counts, {
  '1000GENOMES:phase_3_LWK' => 8,
  '1000GENOMES:phase_3_YRI' => 6,
  '1000GENOMES:phase_3_CLM' => 2,
  '1000GENOMES:phase_3_GWD' => 10,
  '1000GENOMES:phase_3_PUR' => 2,
  '1000GENOMES:phase_3_MSL' => 1,
  '1000GENOMES:phase_3_ASW' => 3,
  '1000GENOMES:phase_3_ESN' => 13,
  '1000GENOMES:phase_3_ACB' => 2,
  '1000GENOMES:phase_3_AFR' => 43,
  '1000GENOMES:phase_3_AMR' => 4
}, 'get_all_population_counts');

is_deeply($h->get_all_population_frequencies, {
  '1000GENOMES:phase_3_LWK' => '0.0404040404040404',
  '1000GENOMES:phase_3_YRI' => '0.0277777777777778',
  '1000GENOMES:phase_3_GWD' => '0.0442477876106195',
  '1000GENOMES:phase_3_CLM' => '0.0106382978723404',
  '1000GENOMES:phase_3_PUR' => '0.00961538461538462',
  '1000GENOMES:phase_3_ASW' => '0.0245901639344262',
  '1000GENOMES:phase_3_MSL' => '0.00588235294117647',
  '1000GENOMES:phase_3_ESN' => '0.0656565656565657',
  '1000GENOMES:phase_3_ACB' => '0.0104166666666667',
  '1000GENOMES:phase_3_AMR' => '0.00576368876080692',
  '1000GENOMES:phase_3_AFR' => '0.0325264750378215'
}, 'get_all_population_frequencies');


$diffs = $h->get_all_diffs();
is_deeply($diffs, [
  {
    'polyphen_prediction' => 'benign',
    'polyphen_score' => '0',
    'diff' => '566R>Q'
  },
  {
    'diff' => '620*>R'
  },
  {
    'polyphen_prediction' => 'benign',
    'polyphen_score' => '0',
    'diff' => '671C>R'
  },
  {
    'polyphen_prediction' => 'probably damaging',
    'polyphen_score' => '0.949',
    'diff' => '678E>A'
  },
  {
    'polyphen_prediction' => 'benign',
    'polyphen_score' => '0',
    'diff' => '819Q>R'
  },
  {
    'polyphen_prediction' => 'probably damaging',
    'polyphen_score' => '0.999',
    'diff' => '899E>L'
  }
], "protein diffs");

($h) = @{$h->get_all_CDSHaplotypes};
$diffs = $h->get_all_diffs();
is(scalar @$diffs, 9, "cds diffs - count");
is($diffs->[0]->{diff}, '1539G>A', "cds diffs - raw name");
is($diffs->[0]->{variation_feature}->variation_name, 'rs114371258', "cds diffs - variation feature");
is($diffs->[-1]->{diff}, '2695GA>TT', "cds diffs - compound");
is($diffs->[-1]->{variation_feature}, undef, "cds diffs - compound no variation feature");


## TranscriptDiplotype tests
my $dts = $c->get_all_CDSDiplotypes;

is(scalar @$dts, 163, "CDS diplotypes - count");

my $dt = (sort {$b->count <=> $a->count} @$dts)[0];

is($dt->name, "ENST00000502692:REF_1697G>A,1858T>C,2011T>C,2456A>G", "CDS diplotype name");
is($dt->count, 636, "CDS diplotype count");
is($dt->frequency, 0.253993610223642, "CDS diplotype freq");
is($dt->container, $c, "CDS diplotype container");
is($dt->transcript, $tr, "CDS diplotype transcript");
is($dt->type, 'cds', "CDS diplotype type");
is_deeply($dt->get_all_population_counts, {
  '1000GENOMES:phase_3_YRI' => 13,
  '1000GENOMES:phase_3_GWD' => 14,
  '1000GENOMES:phase_3_GIH' => 24,
  '1000GENOMES:phase_3_PUR' => 41,
  '1000GENOMES:phase_3_MSL' => 9,
  '1000GENOMES:phase_3_ESN' => 11,
  '1000GENOMES:phase_3_ITU' => 21,
  '1000GENOMES:phase_3_ACB' => 13,
  '1000GENOMES:phase_3_EAS' => 150,
  '1000GENOMES:phase_3_PEL' => 18,
  '1000GENOMES:phase_3_CEU' => 33,
  '1000GENOMES:phase_3_CDX' => 19,
  '1000GENOMES:phase_3_FIN' => 30,
  '1000GENOMES:phase_3_AFR' => 83,
  '1000GENOMES:phase_3_PJL' => 26,
  '1000GENOMES:phase_3_LWK' => 7,
  '1000GENOMES:phase_3_STU' => 29,
  '1000GENOMES:phase_3_IBS' => 35,
  '1000GENOMES:phase_3_CLM' => 32,
  '1000GENOMES:phase_3_EUR' => 164,
  '1000GENOMES:phase_3_SAS' => 124,
  '1000GENOMES:phase_3_ASW' => 16,
  '1000GENOMES:phase_3_CHS' => 35,
  '1000GENOMES:phase_3_TSI' => 31,
  '1000GENOMES:phase_3_MXL' => 24,
  '1000GENOMES:phase_3_GBR' => 35,
  '1000GENOMES:phase_3_CHB' => 34,
  '1000GENOMES:phase_3_BEB' => 24,
  '1000GENOMES:phase_3_AMR' => 115,
  '1000GENOMES:phase_3_KHV' => 23,
  '1000GENOMES:phase_3_JPT' => 39
}, "CDS diplotype pop counts");

is_deeply($dt->get_all_population_frequencies, {
  '1000GENOMES:phase_3_YRI' => '0.12037037037037',
  '1000GENOMES:phase_3_GIH' => '0.233009708737864',
  '1000GENOMES:phase_3_GWD' => '0.123893805309735',
  '1000GENOMES:phase_3_PUR' => '0.394230769230769',
  '1000GENOMES:phase_3_MSL' => '0.105882352941176',
  '1000GENOMES:phase_3_ITU' => '0.205882352941176',
  '1000GENOMES:phase_3_ESN' => '0.111111111111111',
  '1000GENOMES:phase_3_ACB' => '0.135416666666667',
  '1000GENOMES:phase_3_CEU' => '0.333333333333333',
  '1000GENOMES:phase_3_PEL' => '0.211764705882353',
  '1000GENOMES:phase_3_EAS' => '0.297619047619048',
  '1000GENOMES:phase_3_AFR' => '0.125567322239032',
  '1000GENOMES:phase_3_FIN' => '0.303030303030303',
  '1000GENOMES:phase_3_CDX' => '0.204301075268817',
  '1000GENOMES:phase_3_PJL' => '0.270833333333333',
  '1000GENOMES:phase_3_LWK' => '0.0707070707070707',
  '1000GENOMES:phase_3_STU' => '0.284313725490196',
  '1000GENOMES:phase_3_IBS' => '0.327102803738318',
  '1000GENOMES:phase_3_CLM' => '0.340425531914894',
  '1000GENOMES:phase_3_SAS' => '0.253578732106339',
  '1000GENOMES:phase_3_EUR' => '0.326043737574553',
  '1000GENOMES:phase_3_ASW' => '0.262295081967213',
  '1000GENOMES:phase_3_CHS' => '0.333333333333333',
  '1000GENOMES:phase_3_TSI' => '0.289719626168224',
  '1000GENOMES:phase_3_MXL' => '0.375',
  '1000GENOMES:phase_3_GBR' => '0.384615384615385',
  '1000GENOMES:phase_3_AMR' => '0.331412103746398',
  '1000GENOMES:phase_3_BEB' => '0.27906976744186',
  '1000GENOMES:phase_3_CHB' => '0.330097087378641',
  '1000GENOMES:phase_3_JPT' => '0.375',
  '1000GENOMES:phase_3_KHV' => '0.232323232323232'
}, "CDS diplotype pop freqs");


$dts = $c->get_all_ProteinDiplotypes;

is(scalar @$dts, 82, "Protein diplotypes - count");

$dt = (sort {$b->count <=> $a->count} @$dts)[0];

is($dt->name, 'ENSP00000422007:621del{325}_566R>Q,620*>R,671C>R,819Q>R', "Protein diplotype name");
is($dt->count, 676, "Protein diplotype count");
is($dt->frequency, 0.269968051118211, "Protein diplotype freq");
is($dt->container, $c, "Protein diplotype container");
is($dt->transcript, $tr, "Protein diplotype transcript");
is($dt->type, 'protein', "Protein diplotype type");
is_deeply($dt->get_all_population_counts, {
  '1000GENOMES:phase_3_YRI' => 15,
  '1000GENOMES:phase_3_GWD' => 15,
  '1000GENOMES:phase_3_GIH' => 25,
  '1000GENOMES:phase_3_PUR' => 44,
  '1000GENOMES:phase_3_MSL' => 11,
  '1000GENOMES:phase_3_ESN' => 15,
  '1000GENOMES:phase_3_ITU' => 21,
  '1000GENOMES:phase_3_ACB' => 13,
  '1000GENOMES:phase_3_EAS' => 160,
  '1000GENOMES:phase_3_PEL' => 21,
  '1000GENOMES:phase_3_CEU' => 33,
  '1000GENOMES:phase_3_CDX' => 20,
  '1000GENOMES:phase_3_FIN' => 30,
  '1000GENOMES:phase_3_AFR' => 96,
  '1000GENOMES:phase_3_PJL' => 28,
  '1000GENOMES:phase_3_LWK' => 9,
  '1000GENOMES:phase_3_STU' => 30,
  '1000GENOMES:phase_3_IBS' => 35,
  '1000GENOMES:phase_3_CLM' => 34,
  '1000GENOMES:phase_3_EUR' => 167,
  '1000GENOMES:phase_3_SAS' => 129,
  '1000GENOMES:phase_3_ASW' => 18,
  '1000GENOMES:phase_3_CHS' => 38,
  '1000GENOMES:phase_3_TSI' => 33,
  '1000GENOMES:phase_3_MXL' => 25,
  '1000GENOMES:phase_3_GBR' => 36,
  '1000GENOMES:phase_3_CHB' => 37,
  '1000GENOMES:phase_3_BEB' => 25,
  '1000GENOMES:phase_3_AMR' => 124,
  '1000GENOMES:phase_3_KHV' => 24,
  '1000GENOMES:phase_3_JPT' => 41
}, "Protein diplotype pop counts");

is_deeply($dt->get_all_population_frequencies, {
  '1000GENOMES:phase_3_YRI' => '0.138888888888889',
  '1000GENOMES:phase_3_GIH' => '0.242718446601942',
  '1000GENOMES:phase_3_GWD' => '0.132743362831858',
  '1000GENOMES:phase_3_PUR' => '0.423076923076923',
  '1000GENOMES:phase_3_MSL' => '0.129411764705882',
  '1000GENOMES:phase_3_ITU' => '0.205882352941176',
  '1000GENOMES:phase_3_ESN' => '0.151515151515152',
  '1000GENOMES:phase_3_ACB' => '0.135416666666667',
  '1000GENOMES:phase_3_CEU' => '0.333333333333333',
  '1000GENOMES:phase_3_PEL' => '0.247058823529412',
  '1000GENOMES:phase_3_EAS' => '0.317460317460317',
  '1000GENOMES:phase_3_AFR' => '0.145234493192133',
  '1000GENOMES:phase_3_FIN' => '0.303030303030303',
  '1000GENOMES:phase_3_CDX' => '0.21505376344086',
  '1000GENOMES:phase_3_PJL' => '0.291666666666667',
  '1000GENOMES:phase_3_LWK' => '0.0909090909090909',
  '1000GENOMES:phase_3_STU' => '0.294117647058824',
  '1000GENOMES:phase_3_IBS' => '0.327102803738318',
  '1000GENOMES:phase_3_CLM' => '0.361702127659574',
  '1000GENOMES:phase_3_SAS' => '0.263803680981595',
  '1000GENOMES:phase_3_EUR' => '0.332007952286282',
  '1000GENOMES:phase_3_ASW' => '0.295081967213115',
  '1000GENOMES:phase_3_CHS' => '0.361904761904762',
  '1000GENOMES:phase_3_TSI' => '0.308411214953271',
  '1000GENOMES:phase_3_MXL' => '0.390625',
  '1000GENOMES:phase_3_GBR' => '0.395604395604396',
  '1000GENOMES:phase_3_AMR' => '0.357348703170029',
  '1000GENOMES:phase_3_BEB' => '0.290697674418605',
  '1000GENOMES:phase_3_CHB' => '0.359223300970874',
  '1000GENOMES:phase_3_JPT' => '0.394230769230769',
  '1000GENOMES:phase_3_KHV' => '0.242424242424242'
}, "Protein diplotype pop freqs");

done_testing();
