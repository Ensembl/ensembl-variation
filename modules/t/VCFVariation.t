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




###################################
### 1000 genomes with genotypes ###
###################################

# fetch by ID
my $coll = $vca->fetch_by_id('1000genomes_phase1');

# now we need to set the filename_template
my $temp = $coll->filename_template();
$temp =~ s/###t\-root###/$dir/;
$coll->filename_template($temp);

my $slice = $sa->fetch_by_region('chromosome', 2, 45406898, 45406898);
my $vf = $coll->get_all_VariationFeatures_by_Slice($slice)->[0];



## VCFVariationFeature tests
############################

is(ref($vf), 'Bio::EnsEMBL::Variation::VCFVariationFeature', 'VCFVF - check class');

is($vf->variation_name, 'rs11125006', 'VCFVF - variation_name');

is($vf->allele_string, 'T/C', 'VCFVF - allele_string');

is($vf->location_identifier, '2:45406898:T_C:1000genomes', 'VCFVF - location_identifier');

is_deeply($vf->slice, $slice, 'VCFVF - slice');

is_deeply($vf->collection, $coll, 'VCFVF - collection');

is(ref($vf->vcf_record), 'Bio::EnsEMBL::IO::Parser::VCF4Tabix', 'VCFVF - vcf_record class');

is($vf->display_consequence, 'intergenic_variant', 'VCFVF - display_consequence');

is_deeply($vf->$_, [], 'VCFVF - '.$_) for map {'get_all_'.$_.'Variations'} qw(Transcript RegulatoryFeature MotifFeature);



## VCFVariation tests
#####################

my $v = $vf->variation;
is(ref($v), 'Bio::EnsEMBL::Variation::VCFVariation', 'VCFV - check class');

is(ref($v->adaptor), 'Bio::EnsEMBL::Variation::DBSQL::VariationAdaptor', 'VCFV - adaptor class');

is_deeply($v->variation_feature, $vf, 'VCFV - variation_feature');
is_deeply($v->get_all_VariationFeatures, [$vf], 'VCFV - get_all_VariationFeatures');

is($v->name, $vf->variation_name, 'VCFV - name');

is($v->source->name, '1000genomes', 'VCFV - source name');
is($v->source->description, '1000genomes_phase1', 'VCFV - source description');

is_deeply($v->collection, $coll, 'VCFV - collection');

is_deeply($v->vcf_record, $vf->vcf_record, 'VCFV - vcf_record');

is($v->minor_allele, 'T', 'VCFV - minor_allele');
is($v->minor_allele_count, 2, 'VCFV - minor_allele_count');
is(sprintf('%.4g', $v->minor_allele_frequency), 0.3333, 'VCFV - minor_allele_frequency');

is_deeply(
  [ sort map {
    sprintf(
      '%s %s',
      $_->sample->name,
      $_->genotype_string,
    )
  } @{$v->get_all_SampleGenotypes} ],
  [
    '1000GENOMES:phase_1:HG00096 T|C',
    '1000GENOMES:phase_1:HG00097 C|C',
    '1000GENOMES:phase_1:HG00099 C|T'
  ],
  'VCFV - get_all_SampleGenotypes'
);

is_deeply(
  [ sort map {
    sprintf(
      '%s %s %.4g %i',
      $_->population->name,
      $_->allele,
      $_->frequency,
      $_->count
    )
  } @{$v->get_all_Alleles} ],
  [
    '1000GENOMES:phase_1_ALL C 0.6667 4',
    '1000GENOMES:phase_1_ALL T 0.3333 2',
    '1000GENOMES:phase_1_EUR C 0.6667 4',
    '1000GENOMES:phase_1_EUR T 0.3333 2',
    '1000GENOMES:phase_1_GBR C 0.6667 4',
    '1000GENOMES:phase_1_GBR T 0.3333 2'
  ],
  'VCFV - get_all_Alleles'
);

is_deeply(
  [ sort map {
    sprintf(
      '%s %s %.4g %i',
      $_->population->name,
      $_->genotype_string,
      $_->frequency,
      $_->count
    )
  } @{$v->get_all_PopulationGenotypes} ],
  [
    '1000GENOMES:phase_1_ALL C|C 0.3333 1',
    '1000GENOMES:phase_1_ALL C|T 0.6667 2',
    '1000GENOMES:phase_1_EUR C|C 0.3333 1',
    '1000GENOMES:phase_1_EUR C|T 0.6667 2',
    '1000GENOMES:phase_1_GBR C|C 0.3333 1',
    '1000GENOMES:phase_1_GBR C|T 0.6667 2'
  ],
  'VCFV - get_all_PopulationGenotypes'
);




##########################
### ExAC, no genotypes ###
##########################

$coll = $vca->fetch_by_id('ExAC_0.3');

# now we need to set the filename_template
$temp = $coll->filename_template();
$temp =~ s/###t\-root###/$dir/;
$coll->filename_template($temp);
$coll->filename_template =~ /^$dir/;

$slice = $sa->fetch_by_region('chromosome', 11, 66318824, 66318826);
my $exac_vfs = $coll->get_all_VariationFeatures_by_Slice($slice);

# check name generation
is($exac_vfs->[0]->variation_name, '11:66318824:G_C:ExAC', 'ExAC - name generation');

$vf = $exac_vfs->[1];

is(ref($vf), 'Bio::EnsEMBL::Variation::VCFVariationFeature', 'ExAC VCFVF - check class');

is($vf->variation_name, 'rs145769591', 'ExAC VCFVF - variation_name');

is($vf->allele_string, 'G/A/GG', 'ExAC VCFVF - allele_string');

$v = $vf->variation;

is($v->minor_allele, 'A', 'ExAC VCFV - minor_allele');
is($v->minor_allele_count, 11, 'ExAC VCFV - minor_allele_count');
is(sprintf('%.4g', $v->minor_allele_frequency), 9.071e-05, 'ExAC VCFV - minor_allele_frequency');

is_deeply(
  [ sort map {
    sprintf(
      '%s %s %.4g %i',
      $_->population->name,
      $_->allele,
      $_->frequency,
      $_->count
    )
  } @{$v->get_all_Alleles} ],
  [
    'ExAC:AFR A 0.001677 10',
    'ExAC:AFR G 0.9978 5949',
    'ExAC:AFR GG 0.0005032 3',
    'ExAC:ALL A 9.071e-05 11',
    'ExAC:ALL G 0.9999 121248',
    'ExAC:ALL GG 2.474e-05 3',
    'ExAC:AMR A 0 0',
    'ExAC:AMR G 1 4956',
    'ExAC:AMR GG 0 0',
    'ExAC:Adj A 0.0001433 10',
    'ExAC:Adj G 0.9998 69785',
    'ExAC:Adj GG 4.298e-05 3',
    'ExAC:EAS A 0 0',
    'ExAC:EAS G 1 4542',
    'ExAC:EAS GG 0 0',
    'ExAC:FIN A 0 0',
    'ExAC:FIN G 1 3752',
    'ExAC:FIN GG 0 0',
    'ExAC:NFE A 0 0',
    'ExAC:NFE G 1 39012',
    'ExAC:NFE GG 0 0',
    'ExAC:OTH A 0 0',
    'ExAC:OTH G 1 586',
    'ExAC:OTH GG 0 0',
    'ExAC:SAS A 0 0',
    'ExAC:SAS G 1 10988',
    'ExAC:SAS GG 0 0',
  ],
  'ExAC VCFV - get_all_Alleles'
);

is_deeply($v->get_all_SampleGenotypes, [], 'ExAC VCFV - get_all_SampleGenotypes');
is_deeply($v->get_all_PopulationGenotypes, [], 'ExAC VCFV - get_all_PopulationGenotypes');


## consequence stuff

# bypass cons calculation
$slice = $sa->fetch_by_region('chromosome', 11, 66318815, 66318815);
$vf = $coll->get_all_VariationFeatures_by_Slice($slice, 1)->[0];

is($vf->display_consequence, 'intergenic_variant', 'consequence - bypass - display_consequence');
is_deeply($vf->$_, [], 'consequence - bypass - '.$_) for map {'get_all_'.$_.'Variations'} qw(Transcript RegulatoryFeature MotifFeature);

# get for real
$vf = $coll->get_all_VariationFeatures_by_Slice($slice)->[0];
is_deeply([sort @{$vf->consequence_type}], ['3_prime_UTR_variant', 'NMD_transcript_variant', 'synonymous_variant'], 'consequence - consequence_type');

# do a multi fetch for thorough-ness
$slice = $sa->fetch_by_region('chromosome', 11, 66318811, 66318825);
is_deeply(
  [
    map { sprintf('%i %s', $_->seq_region_start, $_->display_consequence) }
    @{$coll->get_all_VariationFeatures_by_Slice($slice)}
  ],
  [
    '66318811 synonymous_variant',
    '66318814 synonymous_variant',
    '66318814 missense_variant',
    '66318815 missense_variant',
    '66318819 missense_variant',
    '66318824 missense_variant'
  ],
  'consequence - multi fetch'
);


## check var class
$slice = $sa->fetch_by_region('chromosome', 11, 66319825, 66319836);

my $vfs = $coll->get_all_VariationFeatures_by_Slice($slice, 1);
is($vfs->[0]->variation->var_class, 'insertion', 'var class - insertion');
is($vfs->[1]->class_SO_term, 'deletion', 'var class - deletion');


## fetch via transcript
# my $tr = $cdb->get_TranscriptAdaptor->fetch_by_stable_id('ENST00000502692');
# my $tva = $vdb->get_TranscriptVariationAdaptor;

# # setup so TVAs get created from our VCF
# $coll->use_as_source(1);
# $vdb->use_vcf(2);

# my $tvas = $tva->fetch_all_by_Transcripts_with_constraint([$tr]);
# $DB::single = 1;

done_testing();
