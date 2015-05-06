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
no warnings 'once';
$Bio::EnsEMBL::Variation::DBSQL::VCFCollectionAdaptor::CONFIG_FILE = $dir.'/th_vcf_config.json';
my $vca = $vdb->get_VCFCollectionAdaptor();
my $vcf_coll = $vca->fetch_all->[0];
my $temp = $vcf_coll->filename_template();
$temp =~ s/###t\-root###/$dir/;
$vcf_coll->filename_template($temp);
$vcf_coll->_add_Populations_to_Individuals();

# get transcript
my $ta = $cdb->get_TranscriptAdaptor;
my $tr = $ta->fetch_by_stable_id('ENST00000502692');

# get genotypes
my @gts;

# we don't want variants in introns
push @gts, @{$vcf_coll->get_all_IndividualGenotypeFeatures_by_Slice($_->feature_Slice, undef, 1)} for @{$tr->get_all_Exons};

my $c = Bio::EnsEMBL::Variation::TranscriptHaplotypeContainer->new($tr, \@gts, $vdb);

my ($i) = grep {$_->name =~ /NA18499/} @{$vcf_coll->get_all_Individuals};


## TESTS
########

my ($h, @h, $counts, $exp, $diffs);

# create container
ok($c->isa('Bio::EnsEMBL::Variation::TranscriptHaplotypeContainer'), "create container");

ok($c->db->isa('Bio::EnsEMBL::Variation::DBSQL::DBAdaptor'), "db method");



# get methods
$h = $c->get_TranscriptHaplotype_by_hex('25f6ab7d0476132b731c630d743d2868');
ok($h->isa('Bio::EnsEMBL::Variation::CDSHaplotype'), "get by hex - isa");
is($h->name, 'ENST00000502692:1697G>A,1818A>G,1858T>C,2011T>C,2033A>C,2456A>G', "get by hex - name");

$h = $c->get_TranscriptHaplotype_by_name('ENST00000502692:1697G>A,1818A>G,1858T>C,2011T>C,2033A>C,2456A>G');
ok($h->isa('Bio::EnsEMBL::Variation::CDSHaplotype'), "get by name - isa");
is($h->hex, '25f6ab7d0476132b731c630d743d2868', "get by name - hex");

@h = @{$c->get_all_TranscriptHaplotypes};
is(scalar @h, 129, "get all TranscriptHaplotypes - count");

@h = @{$c->get_all_CDSHaplotypes};
is(scalar @h, 75, "get all CDSHaplotypes - count");

@h = @{$c->get_all_ProteinHaplotypes};
is(scalar @h, 54, "get all ProteinHaplotypes - count");

@h = @{$c->get_all_TranscriptHaplotypes_by_Individual($i)};
is(scalar @h, 4, "get_all_TranscriptHaplotypes_by_Individual - count");

$DB::single = 1;

@h = @{$c->get_all_CDSHaplotypes_by_Individual($i)};
is(scalar @h, 2, "get_all_CDSHaplotypes_by_Individual - count");

@h = @{$c->get_all_ProteinHaplotypes_by_Individual($i)};
is(scalar @h, 2, "get_all_ProteinHaplotypes_by_Individual - count");

@h = @{$c->get_all_most_frequent_CDSHaplotypes};
ok(scalar @h == 1, "get_all_most_frequent_CDSHaplotypes");
is($h[0]->name, 'ENST00000502692:1697G>A,1858T>C,2011T>C,2456A>G', "get_all_most_frequent_CDSHaplotypes - name");

@h = @{$c->get_all_most_frequent_ProteinHaplotypes};
ok(scalar @h == 1, "get_all_most_frequent_ProteinHaplotypes");
is($h[0]->name, 'ENST00000502692:566R>Q,620*>R,671C>R,819Q>R', "get_all_most_frequent_ProteinHaplotypes - name");



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

@h = @{$c->get_all_ProteinHaplotypes_by_Individual($i)};
$h = $h[1];

is($h->container, $c, 'container');
is($h->type, 'protein', 'type');
is($h->seq, 'MGIGLCVPPARCWKRLTRSSRGQPLERPAATIIIGLSIVERSKQTQRGPVDFSRSHSAPAAEQECKPRISAARSTPAPPSPAVPGFRTTPCQTFTAWCNSHLRKAGTQIENIEEDFRNGLKLMLLLEVISGERLPRPDKGKMRFHKIANVNKALDFIASKGVKLVSIGAEEIVDGNLKMTLGMIWTIILRFAIQDISVEETSAKEGLLLWCQRKTAPYRNVNVQNFHTSWKDGLALCALIHRHRPDLIDYAKLRKDDPIGNLNTAFEVAEKYLDIPKMLDAEDIVNTPKPDEKAIMTYVSCFYHAFAGAEQAETAANRICKVLAVNQENEKLMEEYEKLASELLEWIRRTVPWLENRVGEPSMSAMQRKLEDFRDYRRLHKPPRIQEKCQLEINFNTLQTKLRLSHRPAFMPSEGKLVSDIANAWRGLEQVEKGYEDWLLSEIRRLQRLQHLAEKFRQKASLHEAWTRGKEEMLSQRDYDSALLQEVRALLRRHEAFESDLAAHQDRVEHIAALAQELNELDYHEAASVNSRCQAICDQWDNLGTLTQKRRDALERMEKLLETIDQLQLEFARRAAPFNNWLDGAVEDLQDVWLVHSVEETQSLLTAHDQFKATLPEADRERGAIMGIQGEIQKICQTYGLRPCSTNPYITLSPQDINTKWDMVRKLVPSRDQTLQEALARQQVNERLRRQFAAQANAIGPWIQAKVEEVGRLAAGLAGSLEEQMAGLRQQEQNIINYKTNIDRLEGDHQLLQESLVFDNKHTVYSMEHIRVGWEQLLTSIARTINEVENQVLTRDAKGLSQEQLNEFRASFNHFDRKRNGMMEPDDFRACLISMGYDLGEVEFARIMTMVDPNAAGVVTFQAFIDFMTRETAETDTTEQVVASFKILAGDKNYITPELLRRELPAKQAEYCIRRMVPYKGSGAPAGALDYVAFSSALYGESDL*', 'seq');
is($h->has_indel, 0, 'has_indel');
is($h->hex, 'dff67d1fb7b3812a10cf152d2d9528f7', 'hex');
is_deeply($h->other_hexes, ['fbb87de3331ea4225df6847e8335d622'], 'other_hexes');
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

$DB::single = 1;

done_testing();
