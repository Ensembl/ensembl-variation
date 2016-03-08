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
$vcf_coll->_add_Populations_to_Samples();

# get transcript
my $ta = $cdb->get_TranscriptAdaptor;
my $tr = $ta->fetch_by_stable_id('ENST00000502692');

# get genotypes
my @gts;

# we don't want variants in introns
push @gts, @{$vcf_coll->get_all_SampleGenotypeFeatures_by_Slice($_->feature_Slice, undef, 1)} for @{$tr->get_all_Exons};

my $c = Bio::EnsEMBL::Variation::TranscriptHaplotypeContainer->new(
  -transcript => $tr,
  -genotypes => \@gts,
  # -samples => $vcf_coll->get_all_Samples,
  -db => $vdb,
);

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
is(scalar @h, 115, "get all TranscriptHaplotypes - count");

@h = @{$c->get_all_CDSHaplotypes};
is(scalar @h, 75, "get all CDSHaplotypes - count");

@h = @{$c->get_all_ProteinHaplotypes};
is(scalar @h, 40, "get all ProteinHaplotypes - count");

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
is($c->total_haplotype_count, 4254, "total haplotype count");

$counts = $c->total_population_counts();
$exp = {
  '1000GENOMES:phase_3_YRI' => 216,
  '1000GENOMES:phase_3_GIH' => 160,
  '1000GENOMES:phase_3_GWD' => 224,
  '1000GENOMES:phase_3_PUR' => 174,
  '1000GENOMES:phase_3_MSL' => 170,
  '1000GENOMES:phase_3_ITU' => 162,
  '1000GENOMES:phase_3_ESN' => 196,
  '1000GENOMES:phase_3_ACB' => 190,
  '1000GENOMES:phase_3_PEL' => 82,
  '1000GENOMES:phase_3_CEU' => 162,
  '1000GENOMES:phase_3_EAS' => 840,
  '1000GENOMES:phase_3_FIN' => 172,
  '1000GENOMES:phase_3_CDX' => 130,
  '1000GENOMES:phase_3_AFR' => 1316,
  '1000GENOMES:phase_3_PJL' => 152,
  '1000GENOMES:phase_3_STU' => 148,
  '1000GENOMES:phase_3_LWK' => 198,
  '1000GENOMES:phase_3_IBS' => 180,
  '1000GENOMES:phase_3_CLM' => 154,
  '1000GENOMES:phase_3_EUR' => 842,
  '1000GENOMES:phase_3_SAS' => 748,
  '1000GENOMES:phase_3_ASW' => 122,
  '1000GENOMES:phase_3_CHS' => 184,
  '1000GENOMES:phase_3_TSI' => 184,
  '1000GENOMES:phase_3_MXL' => 98,
  '1000GENOMES:phase_3_GBR' => 144,
  '1000GENOMES:phase_3_CHB' => 176,
  '1000GENOMES:phase_3_AMR' => 508,
  '1000GENOMES:phase_3_BEB' => 126,
  '1000GENOMES:phase_3_KHV' => 174,
  '1000GENOMES:phase_3_JPT' => 176
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
is($h->frequency, 0.0110484250117536, 'frequency');

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
  '1000GENOMES:phase_3_GWD' => '0.0446428571428571',
  '1000GENOMES:phase_3_CLM' => '0.012987012987013',
  '1000GENOMES:phase_3_PUR' => '0.0114942528735632',
  '1000GENOMES:phase_3_ASW' => '0.0245901639344262',
  '1000GENOMES:phase_3_MSL' => '0.00588235294117647',
  '1000GENOMES:phase_3_ESN' => '0.0663265306122449',
  '1000GENOMES:phase_3_ACB' => '0.0105263157894737',
  '1000GENOMES:phase_3_AMR' => '0.0078740157480315',
  '1000GENOMES:phase_3_AFR' => '0.0326747720364742'
}, 'get_all_population_frequencies');

$diffs = $h->get_all_diffs();
is(scalar @$diffs, 6, "protein diffs - count");
is_deeply($diffs->[-1], {
  'diff' => '899E>L',
  'polyphen_prediction' => 'probably damaging',
  'polyphen_score' => 0.999
}, "protein diffs - compound variant");

($h) = @{$h->get_all_CDSHaplotypes};
$diffs = $h->get_all_diffs();
is(scalar @$diffs, 9, "cds diffs - count");
is($diffs->[0]->{diff}, '1539G>A', "cds diffs - raw name");
is($diffs->[0]->{variation_feature}->variation_name, 'rs114371258', "cds diffs - variation feature");
is($diffs->[-1]->{diff}, '2695GA>TT', "cds diffs - compound");
is($diffs->[-1]->{variation_feature}, undef, "cds diffs - compound no variation feature");


## TranscriptDiplotype tests
my $dts = $c->get_all_CDSDiplotypes;

is(scalar @$dts, 162, "CDS diplotypes - count");

my $dt = (sort {$b->count <=> $a->count} @$dts)[0];

is($dt->name, "ENST00000502692:REF_1697G>A,1858T>C,2011T>C,2456A>G", "CDS diplotype name");
is($dt->count, 636, "CDS diplotype count");
is($dt->frequency, 0.29901269393512, "CDS diplotype freq");
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
  '1000GENOMES:phase_3_GIH' => '0.3',
  '1000GENOMES:phase_3_GWD' => '0.125',
  '1000GENOMES:phase_3_PUR' => '0.471264367816092',
  '1000GENOMES:phase_3_MSL' => '0.105882352941176',
  '1000GENOMES:phase_3_ITU' => '0.259259259259259',
  '1000GENOMES:phase_3_ESN' => '0.112244897959184',
  '1000GENOMES:phase_3_ACB' => '0.136842105263158',
  '1000GENOMES:phase_3_CEU' => '0.407407407407407',
  '1000GENOMES:phase_3_PEL' => '0.439024390243902',
  '1000GENOMES:phase_3_EAS' => '0.357142857142857',
  '1000GENOMES:phase_3_AFR' => '0.126139817629179',
  '1000GENOMES:phase_3_FIN' => '0.348837209302326',
  '1000GENOMES:phase_3_CDX' => '0.292307692307692',
  '1000GENOMES:phase_3_PJL' => '0.342105263157895',
  '1000GENOMES:phase_3_LWK' => '0.0707070707070707',
  '1000GENOMES:phase_3_STU' => '0.391891891891892',
  '1000GENOMES:phase_3_IBS' => '0.388888888888889',
  '1000GENOMES:phase_3_CLM' => '0.415584415584416',
  '1000GENOMES:phase_3_SAS' => '0.331550802139037',
  '1000GENOMES:phase_3_EUR' => '0.389548693586698',
  '1000GENOMES:phase_3_ASW' => '0.262295081967213',
  '1000GENOMES:phase_3_CHS' => '0.380434782608696',
  '1000GENOMES:phase_3_TSI' => '0.33695652173913',
  '1000GENOMES:phase_3_MXL' => '0.489795918367347',
  '1000GENOMES:phase_3_GBR' => '0.486111111111111',
  '1000GENOMES:phase_3_AMR' => '0.452755905511811',
  '1000GENOMES:phase_3_BEB' => '0.380952380952381',
  '1000GENOMES:phase_3_CHB' => '0.386363636363636',
  '1000GENOMES:phase_3_JPT' => '0.443181818181818',
  '1000GENOMES:phase_3_KHV' => '0.264367816091954'
}, "CDS diplotype pop freqs");


$dts = $c->get_all_ProteinDiplotypes;

is(scalar @$dts, 81, "Protein diplotypes - count");

$dt = (sort {$b->count <=> $a->count} @$dts)[0];

is($dt->name, 'ENSP00000422007:621del{325}_566R>Q,620*>R,671C>R,819Q>R', "Protein diplotype name");
is($dt->count, 676, "Protein diplotype count");
is($dt->frequency, 0.31781852374236, "Protein diplotype freq");
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
  '1000GENOMES:phase_3_GIH' => '0.3125',
  '1000GENOMES:phase_3_GWD' => '0.133928571428571',
  '1000GENOMES:phase_3_PUR' => '0.505747126436782',
  '1000GENOMES:phase_3_MSL' => '0.129411764705882',
  '1000GENOMES:phase_3_ITU' => '0.259259259259259',
  '1000GENOMES:phase_3_ESN' => '0.153061224489796',
  '1000GENOMES:phase_3_ACB' => '0.136842105263158',
  '1000GENOMES:phase_3_CEU' => '0.407407407407407',
  '1000GENOMES:phase_3_PEL' => '0.51219512195122',
  '1000GENOMES:phase_3_EAS' => '0.380952380952381',
  '1000GENOMES:phase_3_AFR' => '0.145896656534954',
  '1000GENOMES:phase_3_FIN' => '0.348837209302326',
  '1000GENOMES:phase_3_CDX' => '0.307692307692308',
  '1000GENOMES:phase_3_PJL' => '0.368421052631579',
  '1000GENOMES:phase_3_LWK' => '0.0909090909090909',
  '1000GENOMES:phase_3_STU' => '0.405405405405405',
  '1000GENOMES:phase_3_IBS' => '0.388888888888889',
  '1000GENOMES:phase_3_CLM' => '0.441558441558442',
  '1000GENOMES:phase_3_SAS' => '0.344919786096257',
  '1000GENOMES:phase_3_EUR' => '0.39667458432304',
  '1000GENOMES:phase_3_ASW' => '0.295081967213115',
  '1000GENOMES:phase_3_CHS' => '0.41304347826087',
  '1000GENOMES:phase_3_TSI' => '0.358695652173913',
  '1000GENOMES:phase_3_MXL' => '0.510204081632653',
  '1000GENOMES:phase_3_GBR' => '0.5',
  '1000GENOMES:phase_3_AMR' => '0.488188976377953',
  '1000GENOMES:phase_3_BEB' => '0.396825396825397',
  '1000GENOMES:phase_3_CHB' => '0.420454545454545',
  '1000GENOMES:phase_3_JPT' => '0.465909090909091',
  '1000GENOMES:phase_3_KHV' => '0.275862068965517'
}, "Protein diplotype pop freqs");

done_testing();
