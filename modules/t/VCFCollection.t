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
use Test::Exception;
use Bio::EnsEMBL::Test::MultiTestDB;

use FindBin qw($Bin);
use Cwd;

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Variation::VCFCollection;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');

my $vdb = $multi->get_DBAdaptor('variation');
my $cdb = $multi->get_DBAdaptor('core');

my $dir = $multi->curr_dir();

# test constructor
my $c = Bio::EnsEMBL::Variation::VCFCollection->new(
  -id                     => 'test',
  -type                   => 'local',
  -filename_template      => $dir.'/test-genome-DBs/homo_sapiens/variation/test.vcf.gz',
  -chromosomes            => [1,2,3],
  -sample_prefix          => "s_prefix:",
  -population_prefix      => "p_prefix:",
  -sample_populations     => {
    'HG00096' => ['pop1','pop2'],
    'HG00097' => ['pop3'],
    'HG00099' => ['pop4']
},
  -use_seq_region_synonyms => 1,
);


ok($c->id() eq 'test', "id");
ok($c->type() eq 'local', "type");
ok($c->strict_name_match() eq 0, "strict_name_match");
ok($c->filename_template() eq $dir.'/test-genome-DBs/homo_sapiens/variation/test.vcf.gz', "filename_template");
ok($c->sample_prefix() eq "s_prefix:", "sample_prefix");
ok($c->population_prefix() eq "p_prefix:", "population_prefix");
ok($c->use_seq_region_synonyms() eq "1", "use_seq_region_synonyms");
ok($c->tmpdir() eq cwd(), "tmpdir");

# tell it not to use the DB
ok(!$c->use_db(0), "use_db");

# list chromosomes
my $chrs = $c->list_chromosomes();
ok($chrs && scalar @$chrs == 3 && $chrs->[0] eq '1', "list_chromosomes");

# get samples
my $samples = $c->get_all_Samples();
ok($samples && scalar @$samples == 3, "get_all_Samples count 3");
ok($samples->[0]->name eq 's_prefix:HG00096', "get_all_Samples first name is s_prefix:phase_1:HG00096");

# get populations
my $pops = $c->get_all_Populations();

ok($pops && scalar @$pops == 4, "get_all_Populations count 4");
ok($c->has_Population('p_prefix:pop1'), "has_Population p_prefix:pop1");

# test getter/setters
ok(test_getter_setter($c, 'id', 'new_id'), "get/set id");
ok(test_getter_setter($c, 'filename_template', 'new_filename_template'), "get/set filename_template");
ok(test_getter_setter($c, 'sample_prefix', 'new_sample_prefix'), "get/set sample_prefix");
ok(test_getter_setter($c, 'population_prefix', 'new_population_prefix'), "get/set population_prefix");

# test invalid value for type
ok($c->type('remote'), "set type remote");
throws_ok {$c->type('invalid')} qr/Collection type \w+ invalid/, "set invalid type";



# now create a "real" one from config with DB
my $sa = $cdb->get_SliceAdaptor();
ok($sa && $sa->isa('Bio::EnsEMBL::DBSQL::SliceAdaptor'), "get SliceAdaptor");

my $va = $vdb->get_VariationAdaptor;
ok($va && $va->isa('Bio::EnsEMBL::Variation::DBSQL::VariationAdaptor'), "get VariationAdaptor");

ok($vdb->vcf_config_file($dir.'/vcf_config.json') eq $dir.'/vcf_config.json', "DBAdaptor vcf_config_file");
my $vca = $vdb->get_VCFCollectionAdaptor();

ok($vca && $vca->isa('Bio::EnsEMBL::Variation::DBSQL::VCFCollectionAdaptor'), "isa VCFCollectionAdaptor");

# fetch all
my $collections = $vca->fetch_all();

# fetch by ID
my $coll = $vca->fetch_by_id('1000genomes_phase1');
ok($coll && $coll->isa('Bio::EnsEMBL::Variation::VCFCollection'), "fetch_by_id");

# now we need to set the filename_template
my $temp = $coll->filename_template();
$temp =~ s/###t\-root###/$dir/;
$coll->filename_template($temp);
ok($coll->filename_template =~ /^$dir/, "update filename_template");

# get samples
$samples = $coll->get_all_Samples();
ok($samples && scalar @$samples == 3, "get_all_Samples count 3");
ok($samples->[0]->name eq '1000GENOMES:phase_1:HG00096', "get_all_Samples first name is 1000GENOMES:phase_1:HG00096");

# get populations
$pops = $coll->get_all_Populations();
ok($pops && scalar @$pops == 3, "get_all_Populations count 3");
ok($coll->has_Population('1000GENOMES:phase_1_GBR'), "has_Population 1000GENOMES:phase_1_GBR");

# fetch genotypes by VF
my $v  = $va->fetch_by_name('rs7569578');
ok($v && $v->isa('Bio::EnsEMBL::Variation::Variation'), "get variation rs7569578");
my ($vf) = @{$v->get_all_VariationFeatures};
ok($vf && $vf->isa('Bio::EnsEMBL::Variation::VariationFeature'), "get variation feature");

my $gts = $coll->get_all_SampleGenotypeFeatures_by_VariationFeature($vf);
ok($gts && scalar @$gts == 3, "get_all_SampleGenotypeFeatures_by_VariationFeature count 3");
ok($gts->[0]->genotype_string eq 'T|T', "get_all_SampleGenotypeFeatures_by_VariationFeature first genotype T|T");

# fetch genotypes by slice
my $slice = $sa->fetch_by_region('chromosome', 2, 45401130, 45421130);
ok($slice && $slice->isa('Bio::EnsEMBL::Slice'), "get slice");

$gts = $coll->get_all_SampleGenotypeFeatures_by_Slice($slice);
ok($gts && scalar @$gts == 3, "get_all_SampleGenotypeFeatures_by_Slice count 3");
ok($gts->[0]->genotype_string eq 'T|T', "get_all_SampleGenotypeFeatures_by_Slice first genotype T|T");

# fetch LD genotypes by slice
my ($ld_gts, $pos2name) = @{$coll->_get_all_LD_genotypes_by_Slice($slice)};
ok($ld_gts && ref($ld_gts) eq 'HASH', "_get_all_LD_genotypes_by_Slice element 1 is hash");
ok($pos2name && ref($pos2name) eq 'HASH', "_get_all_LD_genotypes_by_Slice element 2 is hash");

ok(scalar keys %$ld_gts == 65, "genotypes hash has 65 position keys");
ok($ld_gts->{45415424} && scalar keys %{$ld_gts->{45415424}} == 3, "genotypes hash pos 45421006 has 3 genotypes");
is_deeply($ld_gts->{45419542}, {
  'HG00096' => 'T|C',
  'HG00099' => 'C|C',
  'HG00097' => 'C|T'
}, "genotypes hash pos 45419542");

is_deeply($pos2name, {
  '45409047' => 'rs111664974',
  '45413788' => 'rs3851322',
  '45404628' => 'rs4953202',
  '45414330' => 'rs11685453',
  '45414280' => 'rs200716478',
  '45409426' => 'rs17392134',
  '45414629' => 'rs12104907',
  '45410572' => 'rs148810371',
  '45404221' => 'rs3914640',
  '45418073' => 'rs62130762',
  '45417945' => 'rs13012034',
  '45402006' => 'rs7591926',
  '45407132' => 'rs4087318',
  '45420335' => 'rs115798714',
  '45419724' => 'rs13025449',
  '45412892' => 'rs13013994',
  '45417972' => 'rs10176710',
  '45406912' => 'rs13425566',
  '45405066' => 'rs3914643',
  '45416527' => 'rs35033699',
  '45411642' => 'rs13392505',
  '45405049' => 'rs13418468',
  '45408893' => 'rs56311136',
  '45401456' => 'rs62133249',
  '45407136' => 'rs116303659',
  '45416594' => 'rs10490353',
  '45407742' => 'rs11903171',
  '45408189' => 'rs4082958',
  '45407137' => 'rs11682694',
  '45410584' => 'rs12613107',
  '45410922' => 'rs7605658',
  '45414283' => 'rs199959267',
  '45402927' => 'rs117686991',
  '45406704' => 'rs35898031',
  '45408678' => 'rs13011860',
  '45414766' => 'rs12105277',
  '45420469' => 'rs150678278',
  '45406940' => 'rs35148156',
  '45418046' => 'rs72799959',
  '45419796' => 'rs66875838',
  '45402571' => 'rs11125005',
  '45419542' => 'rs13025039',
  '45409988' => 'rs17033132',
  '45410449' => 'rs142989295',
  '45402400' => 'rs7595106',
  '45415870' => 'rs60985044',
  '45420082' => 'rs7569588',
  '45406898' => 'rs11125006',
  '45414716' => 'rs72799957',
  '45411443' => 'rs72799955',
  '45411130' => 'rs7569578',
  '45414256' => 'rs199825428',
  '45414221' => 'rs11685398',
  '45408269' => 'rs6720296',
  '45418074' => 'rs13017211',
  '45404640' => 'rs4953203',
  '45415841' => 'rs12612869',
  '45402029' => 'rs3851319',
  '45414377' => 'rs35679497',
  '45420660' => 'rs62130764',
  '45404086' => 'rs11677981',
  '45411188' => 'rs2342702',
  '45414282' => 'rs202033786',
  '45415424' => 'rs17885373',
  '45406008' => 'rs78097421'
}, "pos2name hash");


ok($coll->assembly() eq "GRCh37", "assembly");
ok($coll->source_name() eq "1000genomes", "source name");
ok($coll->source_url() eq "http://www.1000genomes.org", "source URL");

ok($coll->created() eq "1432745640000", "created");
ok($coll->updated() eq "1432745640000", "updated");
ok($coll->is_remapped() eq "1", "is_remapped");

ok($coll->vcf_collection_close, 'close VCF collection filehandle');

# test synonym fetch
$v = $va->fetch_by_name('rs2299222');
ok($v && $v->isa('Bio::EnsEMBL::Variation::Variation'), "get variation rs2299222");
($vf) = @{$v->get_all_VariationFeatures};
ok($vf && $vf->isa('Bio::EnsEMBL::Variation::VariationFeature'), "get variation feature");

$gts = $coll->get_all_SampleGenotypeFeatures_by_VariationFeature($vf);
ok($gts && scalar @$gts == 0, "get_all_SampleGenotypeFeatures_by_VariationFeature without synonyms 0");

$coll->use_seq_region_synonyms(1);
$gts = $coll->get_all_SampleGenotypeFeatures_by_VariationFeature($vf);
ok($gts && scalar @$gts == 3, "get_all_SampleGenotypeFeatures_by_VariationFeature with synonyms 3");

done_testing();
