# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2022] EMBL-European Bioinformatics Institute
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
  -track_name              => "test_track",
  -use_vcf_consequences => 1,
);


ok($c->id() eq 'test', "id");
ok($c->type() eq 'local', "type");
ok($c->strict_name_match() eq 0, "strict_name_match");
ok($c->filename_template() eq $dir.'/test-genome-DBs/homo_sapiens/variation/test.vcf.gz', "filename_template");
ok($c->sample_prefix() eq "s_prefix:", "sample_prefix");
ok($c->population_prefix() eq "p_prefix:", "population_prefix");
ok($c->use_seq_region_synonyms() eq "1", "use_seq_region_synonyms");
ok($c->track_name() eq "test_track", "track_name");
ok($c->use_vcf_consequences() eq "1", "use_vcf_consequences");
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
ok($gts && scalar @$gts == 1119, "get_all_SampleGenotypeFeatures_by_Slice count 1119");
ok($gts->[0]->genotype_string eq 'C|C', "get_all_SampleGenotypeFeatures_by_Slice first genotype C|C");



ok($coll->assembly() eq "GRCh37", "assembly");
ok($coll->source_name() eq "1000genomes", "source name");
ok($coll->source_url() eq "http://www.1000genomes.org", "source URL");


ok($coll->created() eq "1432745640000", "created");
ok($coll->updated() eq "1432745640000", "updated");
ok($coll->is_remapped() eq "1", "is_remapped");
ok($coll->use_vcf_consequences() eq "1", "use vcf consequences");

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

ok($coll->vcf_collection_close, 'close VCF collection filehandle');

## test ESP data query
$v = $va->fetch_by_name('rs80359165');
($vf) = @{$v->get_all_VariationFeatures};
ok($vf && $vf->isa('Bio::EnsEMBL::Variation::VariationFeature'), "get variation feature for rs80359165");

ok($coll->vcf_collection_close, 'close VCF collection filehandle');

$coll = $vca->fetch_by_id('esp_GRCh37');
ok($coll && $coll->isa('Bio::EnsEMBL::Variation::VCFCollection'), "fetch_by_id esp_GRCh37");
$temp = $coll->filename_template();
$temp =~ s/###t\-root###/$dir/;
$coll->filename_template($temp);
$coll->filename_template =~ /^$dir/;
my @alleles = @{$coll->get_all_Alleles_by_VariationFeature($vf)};
is_deeply(
  [
    map {'p:'.$_->population->name.' a:'.$_->allele.' f:'.sprintf("%.4g", $_->frequency).' c:'.$_->count}
    sort {$a->population->name cmp $b->population->name || $a->allele cmp $b->allele}
    @alleles
  ],
  [
    'p:ESP6500:AA a:A f:1 c:4406',
    'p:ESP6500:AA a:C f:0 c:0',
    'p:ESP6500:EA a:A f:0.9999 c:8597',
    'p:ESP6500:EA a:C f:0.0001163 c:1'
  ],
  'get_all_Alleles_by_VariationFeature - freqs and counts ESP rs80359165'
);

my @population_genotypes = @{$coll->get_all_PopulationGenotypes_by_VariationFeature($vf)};

is_deeply(
  [
    map {'p:'.$_->population->name.' genotype:'.$_->genotype_string.' f:'.sprintf("%.4f", $_->frequency).' c:'.$_->count}
    sort {$a->population->name cmp $b->population->name || $a->genotype_string cmp $b->genotype_string}
    @population_genotypes
  ],
  [
    'p:ESP6500:AA genotype:A|A f:1.0000 c:2203',
    'p:ESP6500:EA genotype:A|A f:0.9998 c:4298',
    'p:ESP6500:EA genotype:A|C f:0.0002 c:1'
  ],
  'get_all_PopulationGenotypes_by_VariationFeature - freqs and counts ESP rs80359165'
);

ok($coll->vcf_collection_close, 'close VCF collection filehandle after ESP annotation');

## test get synonyms by chr in offline mode
## chromosomes in VCF file start with chr

$coll = $vca->fetch_by_id('esp_GRCh37_chr_synonym');
$coll->use_db(0);
ok($coll && $coll->isa('Bio::EnsEMBL::Variation::VCFCollection'), "fetch_by_id esp_GRCh37");
$temp = $coll->filename_template();
$temp =~ s/###t\-root###/$dir/;
$coll->filename_template($temp);
$coll->filename_template =~ /^$dir/;
@alleles = @{$coll->get_all_Alleles_by_VariationFeature($vf)};
is_deeply(
  [
    map {'p:'.$_->population->name.' a:'.$_->allele.' f:'.sprintf("%.4g", $_->frequency).' c:'.$_->count}
    sort {$a->population->name cmp $b->population->name || $a->allele cmp $b->allele}
    @alleles
  ],
  [
    'p:ESP6500:AA a:A f:1 c:4406',
    'p:ESP6500:AA a:C f:0 c:0',
    'p:ESP6500:EA a:A f:0.9999 c:8597',
    'p:ESP6500:EA a:C f:0.0001163 c:1'
  ],
  'get_all_Alleles_by_VariationFeature - freqs and counts ESP rs80359165 - test get chr synonym'
);

ok($coll->vcf_collection_close, 'close VCF collection filehandle after ESP annotation');

# Test using AC when AF, AN and AC in the VCF file
# If AC not used and AC calculated by (AF * AN)
# The calculated AC may not match the AC in the VCF due to rounding

# Testing bi-allelic
my $rsid_test = 'rs6720296';
$v = $va->fetch_by_name($rsid_test);
($vf) = @{$v->get_all_VariationFeatures};
ok($vf && $vf->isa('Bio::EnsEMBL::Variation::VariationFeature'), "get variation feature for ${rsid_test}");

$coll = $vca->fetch_by_id('freq_GRCh37');
$coll->use_db(0);
ok($coll && $coll->isa('Bio::EnsEMBL::Variation::VCFCollection'), "fetch_by_id freq_GRCh37");
$temp = $coll->filename_template();
$temp =~ s/###t\-root###/$dir/;
$coll->filename_template($temp);
$coll->filename_template =~ /^$dir/;
@alleles = @{$coll->get_all_Alleles_by_VariationFeature($vf)};

is_deeply(
  [
    map {'p:'.$_->population->name.' a:'.$_->allele.' f:'.sprintf("%.4g", $_->frequency).' c:'.$_->count}
    sort {$a->population->name cmp $b->population->name || $a->allele cmp $b->allele}
    @alleles
  ],
  [
    'p:test_freq a:A f:0.6251 c:9489',
    'p:test_freq a:C f:0.3749 c:5691',
  ],
  "get_all_Alleles_by_VariationFeature - freqs and counts test_freq $rsid_test - test rounding"
);

# Testing multi-allelic
$rsid_test='rs370058043';
$v = $va->fetch_by_name($rsid_test);
($vf) = @{$v->get_all_VariationFeatures};
ok($vf && $vf->isa('Bio::EnsEMBL::Variation::VariationFeature'), "get variation feature for ${rsid_test}");

@alleles = @{$coll->get_all_Alleles_by_VariationFeature($vf)};

is_deeply(
  [
    map {'p:'.$_->population->name.' a:'.$_->allele.' f:'.sprintf("%.4g", $_->frequency).' c:'.$_->count}
    sort {$a->population->name cmp $b->population->name || $a->allele cmp $b->allele}
    @alleles
  ],
  [
    'p:test_freq a:C f:8.288e-05 c:6',
    'p:test_freq a:G f:0.9998 c:72376',
    'p:test_freq a:T f:0.0001381 c:10',
  ],
  "get_all_Alleles_by_VariationFeature - freqs and counts test_freq $rsid_test - test rounding"
);

ok($coll->vcf_collection_close, 'close VCF collection filehandle after freq_test');


## test exac info stuff

# fetch by ID
$coll = $vca->fetch_by_id('ExAC_0.3');
ok($coll && $coll->isa('Bio::EnsEMBL::Variation::VCFCollection'), "fetch_by_id");

# now we need to set the filename_template
$temp = $coll->filename_template();
$temp =~ s/###t\-root###/$dir/;
$coll->filename_template($temp);
$coll->filename_template =~ /^$dir/;

# check list_chromosomes works with no "chromosomes" field in JSON
delete $coll->{chromosomes};
is_deeply(
  $coll->list_chromosomes,
  ["11"],
  'list_chromosomes repopulates from tabix index'
);

($vf) = @{$va->fetch_by_name('rs192076014')->get_all_VariationFeatures};
@alleles = @{$coll->get_all_Alleles_by_VariationFeature($vf)};

is(scalar @alleles, 18, 'get_all_Alleles_by_VariationFeature - count');
is_deeply(
  [
    map {'p:'.$_->population->name.' a:'.$_->allele.' f:'.sprintf("%.4f", $_->frequency).' c:'.$_->count}
    sort {$a->population->name cmp $b->population->name || $a->allele cmp $b->allele}
    @alleles
  ],
  [
    'p:ExAC:AFR a:A f:0.0034 c:22',
    'p:ExAC:AFR a:G f:0.9966 c:6436',
    'p:ExAC:ALL a:A f:0.0002 c:24',
    'p:ExAC:ALL a:G f:0.9998 c:121312',
    'p:ExAC:AMR a:A f:0.0002 c:1',
    'p:ExAC:AMR a:G f:0.9998 c:5837',
    'p:ExAC:Adj a:A f:0.0003 c:24',
    'p:ExAC:Adj a:G f:0.9997 c:76724',
    'p:ExAC:EAS a:A f:0.0000 c:0',
    'p:ExAC:EAS a:G f:1.0000 c:5244',
    'p:ExAC:FIN a:A f:0.0000 c:0',
    'p:ExAC:FIN a:G f:1.0000 c:4094',
    'p:ExAC:NFE a:A f:0.0000 c:0',
    'p:ExAC:NFE a:G f:1.0000 c:42906',
    'p:ExAC:OTH a:A f:0.0016 c:1',
    'p:ExAC:OTH a:G f:0.9984 c:619',
    'p:ExAC:SAS a:A f:0.0000 c:0',
    'p:ExAC:SAS a:G f:1.0000 c:11588'
  ],
  'get_all_Alleles_by_VariationFeature - freqs and counts'
);

# test with a VCF file which contains a duplicated rows for rs192076014
$coll = $vca->fetch_by_id('ExAC_0.3_with_duplicated_row');
ok($coll && $coll->isa('Bio::EnsEMBL::Variation::VCFCollection'), "fetch_by_id");

# now we need to set the filename_template
$temp = $coll->filename_template();
$temp =~ s/###t\-root###/$dir/;
$coll->filename_template($temp);
$coll->filename_template =~ /^$dir/;

@alleles = @{$coll->get_all_Alleles_by_VariationFeature($vf)};

is(scalar @alleles, 18, 'get_all_Alleles_by_VariationFeature - count from VCF with duplicated row');
is_deeply(
  [
    map {'p:'.$_->population->name.' a:'.$_->allele.' f:'.sprintf("%.4f", $_->frequency).' c:'.$_->count}
    sort {$a->population->name cmp $b->population->name || $a->allele cmp $b->allele}
    @alleles
  ],
  [
    'p:ExAC:AFR a:A f:0.0034 c:22',
    'p:ExAC:AFR a:G f:0.9966 c:6436',
    'p:ExAC:ALL a:A f:0.0002 c:24',
    'p:ExAC:ALL a:G f:0.9998 c:121312',
    'p:ExAC:AMR a:A f:0.0002 c:1',
    'p:ExAC:AMR a:G f:0.9998 c:5837',
    'p:ExAC:Adj a:A f:0.0003 c:24',
    'p:ExAC:Adj a:G f:0.9997 c:76724',
    'p:ExAC:EAS a:A f:0.0000 c:0',
    'p:ExAC:EAS a:G f:1.0000 c:5244',
    'p:ExAC:FIN a:A f:0.0000 c:0',
    'p:ExAC:FIN a:G f:1.0000 c:4094',
    'p:ExAC:NFE a:A f:0.0000 c:0',
    'p:ExAC:NFE a:G f:1.0000 c:42906',
    'p:ExAC:OTH a:A f:0.0016 c:1',
    'p:ExAC:OTH a:G f:0.9984 c:619',
    'p:ExAC:SAS a:A f:0.0000 c:0',
    'p:ExAC:SAS a:G f:1.0000 c:11588'
  ],
  'get_all_Alleles_by_VariationFeature - freqs and counts from VCF with duplicated row'
);


# back to VCF file without duplicated rows 
$coll = $vca->fetch_by_id('ExAC_0.3');
ok($coll && $coll->isa('Bio::EnsEMBL::Variation::VCFCollection'), "fetch_by_id");

# now we need to set the filename_template
$temp = $coll->filename_template();
$temp =~ s/###t\-root###/$dir/;
$coll->filename_template($temp);
$coll->filename_template =~ /^$dir/;

# test one which has an allele not found in dbSNP entry
($vf) = @{$va->fetch_by_name('rs145769591')->get_all_VariationFeatures};
is_deeply(
  [
    map {keys %{$_->{_missing_alleles}}}
    grep {$_->population->name eq 'ExAC:AFR'}
    @{$coll->get_all_Alleles_by_VariationFeature($vf)}
  ],
  ['GG', 'GG'],
  'get_all_Alleles_by_VariationFeature - _missing_alleles'
);

# test one that merges across VCF lines
($vf) = @{$va->fetch_by_name('rs547901734')->get_all_VariationFeatures};
@alleles = @{$coll->get_all_Alleles_by_VariationFeature($vf)};
is_deeply(
  [
    map {'p:'.$_->population->name.' a:'.$_->allele.' f:'.sprintf("%.4f", $_->frequency).' c:'.$_->count}
    sort {$a->population->name cmp $b->population->name || $a->allele cmp $b->allele}
    @alleles
  ],
  [
    'p:ExAC:AFR a:A f:0.0000 c:0',
    'p:ExAC:AFR a:C f:0.0003 c:2',
    'p:ExAC:AFR a:G f:0.9997 c:6220',
    'p:ExAC:ALL a:A f:0.0001 c:12',
    'p:ExAC:ALL a:C f:0.0000 c:2',
    'p:ExAC:ALL a:G f:0.9999 c:121286',
    'p:ExAC:AMR a:A f:0.0004 c:2',
    'p:ExAC:AMR a:C f:0.0000 c:0',
    'p:ExAC:AMR a:G f:0.9996 c:5416',
    'p:ExAC:Adj a:A f:0.0001 c:11',
    'p:ExAC:Adj a:C f:0.0000 c:2',
    'p:ExAC:Adj a:G f:0.9998 c:73547',
    'p:ExAC:EAS a:A f:0.0000 c:0',
    'p:ExAC:EAS a:C f:0.0000 c:0',
    'p:ExAC:EAS a:G f:1.0000 c:4964',
    'p:ExAC:FIN a:A f:0.0000 c:0',
    'p:ExAC:FIN a:C f:0.0000 c:0',
    'p:ExAC:FIN a:G f:1.0000 c:3954',
    'p:ExAC:NFE a:A f:0.0002 c:9',
    'p:ExAC:NFE a:C f:0.0000 c:0',
    'p:ExAC:NFE a:G f:0.9998 c:41113',
    'p:ExAC:OTH a:A f:0.0000 c:0',
    'p:ExAC:OTH a:C f:0.0000 c:0',
    'p:ExAC:OTH a:G f:1.0000 c:600',
    'p:ExAC:SAS a:A f:0.0000 c:0',
    'p:ExAC:SAS a:C f:0.0000 c:0',
    'p:ExAC:SAS a:G f:1.0000 c:11280'
  ],
  'get_all_Alleles_by_VariationFeature - merge across VCF lines'
);


## test dbsnp which has ref freq included
($vf) = @{$va->fetch_by_name('rs192076014')->get_all_VariationFeatures};
$coll = $vca->fetch_by_id('dbsnp');
ok($coll && $coll->isa('Bio::EnsEMBL::Variation::VCFCollection'), "fetch_by_id");

# now we need to set the filename_template
$temp = $coll->filename_template();
$temp =~ s/###t\-root###/$dir/;
$coll->filename_template($temp);
$coll->filename_template =~ /^$dir/;

is_deeply(
  [map {'a:'.$_->allele.' f:'.$_->frequency} sort {$a->allele cmp $b->allele} @{$coll->get_all_Alleles_by_VariationFeature($vf)}],
  [
    'a:A f:0.00120201',
    'a:G f:0.998798'
  ],
  'get_all_Alleles_by_VariationFeature - dbSNP uses ref_freq_index()'
);

# Test get_all_clinical_significance_states() for VCF files
$coll = $vca->fetch_by_id('ClnSig');
ok($coll && $coll->isa('Bio::EnsEMBL::Variation::VCFCollection'), "fetch_by_id");

my $temp = $coll->filename_template();
$temp =~ s/###t\-root###/$dir/;
$coll->filename_template($temp);
ok($coll->filename_template =~ /^$dir/, "update filename_template");

$slice = $sa->fetch_by_region('toplevel', 1, 10, 20);
my $dont_fetch_vf_overlaps=1;
my $vfs = $coll->get_all_VariationFeatures_by_Slice($slice, $dont_fetch_vf_overlaps);

ok($vfs->[0]->get_all_clinical_significance_states()->[0] eq 'likely benign', 'get_all_clinical_significance_states - obtain single clinical significance entry');
ok(scalar (@{$vfs->[1]->get_all_clinical_significance_states()}) eq 1, 'get_all_clinical_significance_states - ignore upsupported clinical significance entry');
ok(scalar (@{$vfs->[2]->get_all_clinical_significance_states()}) eq 2, 'get_all_clinical_significance_states - obtain multiple clinical significance entries');

done_testing();
