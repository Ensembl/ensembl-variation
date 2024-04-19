# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2024] EMBL-European Bioinformatics Institute
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
use Bio::EnsEMBL::Variation::Source;
use Bio::EnsEMBL::Test::MultiTestDB;

use_ok('Bio::EnsEMBL::Variation::Variation');
use_ok('Bio::EnsEMBL::Variation::VariationFeature');


my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');
my $vdb = $multi->get_DBAdaptor('variation');
my $cdb = $multi->get_DBAdaptor('core');

my $vfa = $vdb->get_VariationFeatureAdaptor();
my $va  = $vdb->get_VariationAdaptor();
my $sa  = $cdb->get_SliceAdaptor();


# test constructor

## need source object

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


my $chr   = 18;
my $start = 100000;
my $end   = 100000;
my $strand = 1;
my $vname = $v->name();
my $map_weight = 1;
my $allele_str = 'A/T';
my $ancestral_allele = 'A';
my $is_somatic = 0;
my $minor_allele = 'A';
my $minor_allele_frequency = 0.1;
my $minor_allele_count = 10;
my $clin_sign = ['likely pathogenic','pathogenic'];
my $consequence = Bio::EnsEMBL::Variation::OverlapConsequence->new(-SO_term => 'stop_gained');

my $slice = $sa->fetch_by_region('chromosome',$chr,$start,$end,$strand);

my $vf = Bio::EnsEMBL::Variation::VariationFeature->new
  (-seq_region_name => $chr,
   -start => $start,
   -end   => $end,
   -slice => $slice,
   -strand => $strand,
   -variation_name => $vname,
   -map_weight => $map_weight,
   -allele_string => $allele_str,
   -ancestral_allele => $ancestral_allele,
   -variation => $v,
   -source => $source,
   -is_somatic => $is_somatic,
   -minor_allele  => $minor_allele ,
   -minor_allele_frequency  => $minor_allele_frequency ,
   -minor_allele_count => $minor_allele_count,
   -clinical_significance => $clin_sign,
   -overlap_consequences => [$consequence],
);

ok($vf->seq_region_name() == $chr,      "get chromosome name");
ok($vf->start()   == $start,            "get start");
ok($vf->end()     == $end,              "get end");
ok($vf->strand()  == $strand,           "get strand");
ok($vf->variation_name() eq $vname,     "get variation_name");
ok($vf->name() eq $vname,               "get name");
ok($vf->map_weight() == $map_weight,    "get map_weight");
ok($vf->allele_string() eq $allele_str, "get allele");
ok($vf->ancestral_allele() eq $ancestral_allele, "get ancestral_allele");
ok($vf->ref_allele_string() eq 'A',     "get ref allele");
is_deeply($vf->alt_alleles, ['T'],      "get alt alleles");
ok($vf->display_id() eq $vname,         "display_name");
ok($vf->source_name() eq $source_name,  "source");
ok($vf->source_version() eq $source_version, "source version");
ok($vf->length() == 1,                  "length");
ok($vf->is_somatic() eq $is_somatic,    "is_somatic");
ok($vf->minor_allele() eq $minor_allele, "minor allele");
ok($vf->minor_allele_frequency() == $minor_allele_frequency , "minor allele freq");
ok($vf->minor_allele_count() == $minor_allele_count,  "minor allele count");
ok($vf->get_Variation_dbID() == 12345,  "get variation db id");
ok($vf->is_reference,                   "check that the VF's slice is a reference");
ok($vf->consequence_type->[0] eq $consequence->SO_term, "get the consequence");

ok($vf->add_evidence_value("Cited"), 'add a permitted evidence value');
my $oc = Bio::EnsEMBL::Variation::OverlapConsequence->new(-SO_term => 'missense_variant');
ok($vf->add_OverlapConsequence($oc), 'add_OverlapConsequence');
ok($vf->get_all_evidence_values()->[0] eq 'Cited', 'get_all_evidence_values');
ok($vf->get_all_clinical_significance_states()->[0] eq $clin_sign->[0], 'get_all_clinical_significance_states');
is($vf->feature_so_acc, 'SO:0001060', 'VariationFeature feature SO acc is correct (sequence variant)');
is($vf->feature_so_term, 'sequence_variant', 'VariationFeature feature SO term is correct (sequence variant)');

# test getter/setters

my $v2 = Bio::EnsEMBL::Variation::Variation->new(-name => 'rs12311',
                                                 -source => $source);

ok(test_getter_setter($vf, 'variation', $v2),            'set new variation object');
ok(test_getter_setter($vf, 'map_weight', 4),             'set new map_weight');
ok(test_getter_setter($vf, 'allele_string', 'T/G'),    'set new allele_string');
ok(test_getter_setter($vf, 'ancestral_allele', 'T'),    'set new ancestral_allele');
ok(test_getter_setter($vf, 'variation_name', 'rs21431'), 'set new variation name');
ok(test_getter_setter($vf, 'flank_match', '1'),          'set new flank_match');


my $bak = $vf->allele_string;
$vf->allele_string('T/G/C');
ok($vf->ref_allele_string() eq 'T',     "get ref allele after change");
is_deeply($vf->alt_alleles, ['G', 'C'], "get alt alleles after change");
$vf->allele_string($bak);


# test ambiguity code
ok($vf->ambig_code eq 'W', "ambiguity code");

# test variation class
ok($vf->var_class eq 'SNP', "var class");



## VariationFeature_to_VCF_record
#################################

my $sr_start = $vf->seq_region_start;

is_deeply(
  $vf->to_VCF_record(),
  [$chr, $sr_start, $vname, 'A', 'T', '.', '.', "END=$end"],
  'to_VCF_record'
);

$vf->strand(-1);
is_deeply(
  $vf->to_VCF_record(),
  [$chr, $sr_start, $vname, 'T', 'A', '.', '.', "END=$end"],
  'to_VCF_record - rev strand'
);
$vf->strand($strand);

$vf->allele_string('A/G/T');
is_deeply(
  $vf->to_VCF_record(),
  [$chr, $sr_start, $vname, 'A', 'G,T', '.', '.', "END=$end"],
  'to_VCF_record - multiple alts'
);

$vf->allele_string('AG/CT');
is_deeply(
  $vf->to_VCF_record(),
  [$chr, $sr_start, $vname, 'AG', 'CT', '.', '.', "END=$end"],
  'to_VCF_record - balanced non-SNP'
);

$vf->allele_string('A/-');
is_deeply(
  $vf->to_VCF_record(),
  [$chr, $sr_start - 1, $vname, 'NA', 'N', '.', '.', "END=$end"],
  'to_VCF_record - deletion'
);

$vf->allele_string('-/A');
is_deeply(
  $vf->to_VCF_record(),
  [$chr, $sr_start - 1, $vname, 'N', 'NA', '.', '.', "END=$end"],
  'to_VCF_record - insertion'
);

$vf->allele_string('A/-/G');
is_deeply(
  $vf->to_VCF_record(),
  [$chr, $sr_start - 1, $vname, 'NA', 'N,NG', '.', '.', "END=$end"],
  'to_VCF_record - mixed'
);


$vf->allele_string('HGMD_MUTATION');
$vf->{class_SO_term} = 'SNV';
is_deeply(
  $vf->to_VCF_record(),
  [$chr, $sr_start, $vname, 'N', 'N', '.', '.', "END=$end"],
  'to_VCF_record - unknown alleles SNV'
);

$vf->{class_SO_term} = 'insertion';
is_deeply(
  $vf->to_VCF_record(),
  [$chr, $sr_start - 1, $vname, 'N', '<INS>', '.', '.', "END=$end"],
  'to_VCF_record - unknown alleles insertion'
);

$vf->{class_SO_term} = 'deletion';
is_deeply(
  $vf->to_VCF_record(),
  [$chr, $sr_start - 1, $vname, 'NN', 'N', '.', '.', "END=$end"],
  'to_VCF_record - unknown alleles deletion'
);

$vf->{class_SO_term} = 'sequence_alteration';
is_deeply(
  $vf->to_VCF_record(),
  [],
  'to_VCF_record - unknown alleles sequence_alteration'
);

my $fully_justified_allele_str = 'ACGTGGACG/ACG/ACGTGGACGTGGACG';
$vf->allele_string($fully_justified_allele_str);
is_deeply(
  $vf->to_VCF_record(),
  [$chr, $sr_start, $vname, 'ACGTGGA', 'A,ACGTGGACGTGGA', '.', '.', "END=$end"],
  'to_VCF_record - fully justified allele string clipped'
);


is_deeply(
  $vf->to_VCF_record(1),
  [$chr, $sr_start, $vname, 'ACGTGGACG', 'ACG,ACGTGGACGTGGACG', '.', '.', "END=$end"],
  'to_VCF_record - fully justified allele string  not clipped'
);


$vf->allele_string($allele_str);


# test overridden methods
my $sr_coord = $start + $start - 1;
is($vf->seq_region_start, $sr_coord, 'seq_region_start');
is($vf->{seq_region_start}, $sr_coord, 'seq_region_start cached');
is($vf->seq_region_start(10), 10, 'seq_region_start set');

is($vf->seq_region_end, $sr_coord, 'seq_region_end');
is($vf->{seq_region_end}, $sr_coord, 'seq_region_end cached');
is($vf->seq_region_end(10), 10, 'seq_region_end set');

my $vf_transformed = $vf->transform('contig');
ok(
  (
    !exists($vf->{seq_region_start}) &&
    !exists($vf->{seq_region_end}) &&
    !exists($vf_transformed->{seq_region_start}) &&
    !exists($vf_transformed->{seq_region_end})
  ),
  'transform deletes seq_region_start and seq_region_end'
);

# repopulate cache
$vf->seq_region_start;
$vf->seq_region_end;
is($vf->{seq_region_start}, $sr_coord, 'seq_region_start re-cached');
is($vf->{seq_region_end}, $sr_coord, 'seq_region_end re-cached');

my $exp_slice = $slice->expand(10, 10);
my $vf_transferred = $vf->transfer($exp_slice);
ok(
  (
    !exists($vf->{seq_region_start}) &&
    !exists($vf->{seq_region_end}) &&
    !exists($vf_transferred->{seq_region_start}) &&
    !exists($vf_transferred->{seq_region_end})
  ),
  'transfer deletes seq_region_start and seq_region_end'
);
is($vf_transferred->seq_region_start, $vf->seq_region_start, 'seq_region_start after transfer');
is($vf_transferred->seq_region_end, $vf->seq_region_end, 'seq_region_end after transfer');

# test convert to SNP
#ok($vf->convert_to_SNP, 'convert to SNP'); # Need the ensembl-external repository


# test get all VariationSets
my $var3 = $va->fetch_by_name('rs80359159');
my $vf3 = $var3->get_all_VariationFeatures()->[0];
my $vss = $vf3->get_all_VariationSets();
my $vs_name = 'clinically associated';
ok( scalar (grep { $_->name eq $vs_name } @$vss) == 1, 'get all VariationSets');

# get all Alleles
my $als = $vf3->get_all_Alleles();
ok(scalar(@$als) == 7 && $als->[0]->allele eq 'C', 'get all Alleles');

# test get all PopulationGenotypes
my $var4 = $va->fetch_by_name('rs2255888');
my $vf4 = $var4->get_all_VariationFeatures()->[0];
my $pgs = $vf4->get_all_PopulationGenotypes();
ok($pgs->[0]->genotype_string eq 'C|C' && $pgs->[0]->count == 98, 'get all PopulationGenotypes');

# test get all sources
ok($vf4->get_all_sources->[0] eq 'dbSNP', 'get_all_sources');

# test get all RegulatoryFeatureVariations
my $var5 = $va->fetch_by_name('rs187207343');
my $vf5 = $var5->get_all_VariationFeatures()->[0];
my $rfs = $vf5->get_all_RegulatoryFeatureVariations;
ok($rfs->[0]->regulatory_feature_stable_id eq 'ENSR00000000637', 'get_all_RegulatoryFeatureVariations');
my $regulatory_feature = $rfs->[0]->regulatory_feature;
$rfs = $vf5->get_all_RegulatoryFeatureVariations([$regulatory_feature]);
ok($rfs->[0]->regulatory_feature_stable_id eq 'ENSR00000000637', 'get_all_RegulatoryFeatureVariations, regulatory_feature');

# test get all MotifFeatureVariations
my $var6 = $va->fetch_by_name('rs182313188');
my $vf6 = $var6->get_all_VariationFeatures()->[0];
my $mfvs = $vf6->get_all_MotifFeatureVariations;
ok($mfvs->[0]->feature_stable_id eq 'ENSR00000636355', 'get_all_MotifFeatureVariations');
my $motif_feature = $mfvs->[0]->motif_feature;
$mfvs = $vf6->get_all_MotifFeatureVariations([$motif_feature]);
ok($mfvs->[0]->feature_stable_id eq 'ENSR00000636355', 'get_all_MotifFeatureVariations, motif_feature');

# test spdi genomic
# To do: SPDI annotation should be updated from chrom number to synonym, i.e: 13 -> NC_000013.10
my $spdi_notation_3 = $vf3->spdi_genomic();
ok($spdi_notation_3->{'T'} eq 'NC_000013.10:32954008:C:T', 'SPDI genomic valid substitution - alt allele 1');
ok($spdi_notation_3->{'G'} eq 'NC_000013.10:32954008:C:G', 'SPDI genomic valid substitution - alt allele 2');
my $var7 = $va->fetch_by_name('rs35370278');
my $vf7_spdi = $var7->get_all_VariationFeatures()->[0];
my $spdi_notation_7 = $vf7_spdi->spdi_genomic();
ok($spdi_notation_7->{'-'} eq 'NC_000011.9:66315558:G:', 'SPDI genomic valid deletion');
my $spdi_notation_7_ref_allele = $vf7_spdi->spdi_genomic(1);
ok($spdi_notation_7_ref_allele->{'G'} eq 'NC_000011.9:66315558:G:G', 'SPDI genomic valid deletion - ref allele');
throws_ok {$vf7_spdi->spdi_genomic(2); } qr/Include reference allele must be a numeric value '1' or '0'./, 'Throw invalid input'; 
my $var8 = $va->fetch_by_name('rs370045702');
my $vf8_spdi = $var8->get_all_VariationFeatures()->[0];
my $spdi_notation_8 = $vf8_spdi->spdi_genomic();
ok($spdi_notation_8->{'AA'} eq 'NC_000011.9:66317203::AA', 'SPDI genomic valid insertion');
my $var9 = $va->fetch_by_name('rs35794957');
my $vf9_spdi = $var9->get_all_VariationFeatures()->[0];
my $spdi_notation_9 = $vf9_spdi->spdi_genomic();
ok($spdi_notation_9->{'CA'} eq 'NC_000011.9:66321302:TG:CA', 'SPDI genomic valid indel');
my $var10 = $va->fetch_by_name('rs1480172069');
my $vf10_spdi = $var10->get_all_VariationFeatures()->[0];
my $spdi_notation_10 = $vf10_spdi->spdi_genomic();
ok($spdi_notation_10->{'-'} eq 'NC_000019.9:48835629:ATTT:', 'SPDI genomic valid deletion - reverse strand'); 
my $var11 = $va->fetch_by_name('rs1321600644');
my $vf11_spdi = $var11->get_all_VariationFeatures()->[0];
my $spdi_notation_11 = $vf11_spdi->spdi_genomic();
ok($spdi_notation_11->{'-'} eq 'NC_000012.11:101997654:25:', 'SPDI genomic valid deletion >20 bp'); 
my $var12 = $va->fetch_by_name('rs8192742');
my $vf12_spdi = $var12->get_all_VariationFeatures()->[0];
my $spdi_notation_12 = $vf12_spdi->spdi_genomic();
ok($spdi_notation_12->{'TTTTTTTTTT'} eq 'NC_000013.10:76134860:TTTTTTTTT:TTTTTTTTTT', 'SPDI tandem repeat');

#test deprecated methods
print "\n## Test deprecated methods ##\n";

# LD data
# VCF
my $dir = $multi->curr_dir();
ok($vdb->vcf_config_file($dir.'/ld_vcf_config.json') eq $dir.'/ld_vcf_config.json', "DBAdaptor vcf_config_file");
my $vca = $vdb->get_VCFCollectionAdaptor();
my $coll = $vca->fetch_by_id('ld');
# now we need to set the filename_template
my $temp = $coll->filename_template();
$temp =~ s/###t\-root###/$dir/;
$coll->filename_template($temp);
$vfa->db->use_vcf(1);
my $vf7 = $vfa->fetch_by_dbID(3854101);
my @LD_populations = @{$vf7->get_all_LD_Populations};
ok(scalar @LD_populations == 2, 'get_all_LD_Populations use_vcf = 1');
$vfa->db->use_vcf(2);
@LD_populations = @{$vf7->get_all_LD_Populations};
ok(scalar @LD_populations == 1, 'get_all_LD_Populations use_vcf = 2');
$vfa->db->use_vcf(0);
@LD_populations = @{$vf7->get_all_LD_Populations};
ok(scalar @LD_populations == 1, 'get_all_LD_Populations use_vcf = 0');


# get_all_highest_frequency_minor_Alleles
$vfa->db->use_vcf(1);
my $hpmaf_alleles = $vf7->get_all_highest_frequency_minor_Alleles();
is(scalar @$hpmaf_alleles, 2, 'get_all_highest_frequency_minor_Alleles - count');
is(sprintf("%.4f", $hpmaf_alleles->[0]->frequency), 0.1148, 'get_all_highest_frequency_minor_Alleles - frequency');

is_deeply(
  [map {$_->allele} @$hpmaf_alleles],
  [qw(C C)],
  'get_all_highest_frequency_minor_Alleles - allele'
);
is_deeply(
  [sort map {$_->population->name} @$hpmaf_alleles],
  [qw(1000GENOMES:phase_1_AFR 1000GENOMES:phase_1_ASW)],
  'get_all_highest_frequency_minor_Alleles - population'
);

$vfa->db->use_vcf(0);


# test fake variation feature adapter and fake variation feature
{
  my $vfa = Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor->new_fake('homo_sapiens');
  my $tmp_vf_fs = Bio::EnsEMBL::Variation::VariationFeature->new_fast({
    start          => 100000,
    end            => 100005,
    allele_string  => 'A/-',
    strand         => 1,
    map_weight     => 1,
    adaptor        => $vfa,
    slice          => $slice,
    seq_region_start => $slice->seq_region_start,
    seq_region_end   => $slice->seq_region_end
  });
  my $conseq = join ",", @{$tmp_vf_fs->consequence_type};
  ok($conseq eq 'intergenic_variant', "fake VariationFeatureAdaptor and VariationFeature");
}

# test most_severe_OverlapConsequence for consequences with the same rank
my $oc_1 = Bio::EnsEMBL::Variation::OverlapConsequence->new(
             -SO_term => 'splice_donor_variant',
             -rank    => 3);

my $oc_2 = Bio::EnsEMBL::Variation::OverlapConsequence->new(
              -SO_term => 'splice_acceptor_variant',
              -rank    => 3);

my $vf_msc = Bio::EnsEMBL::Variation::VariationFeature->new
  ( -overlap_consequences => [$oc_1, $oc_2]);

my $msc_2_expected = 'splice_acceptor_variant';
my $msc_2 = $vf_msc->most_severe_OverlapConsequence();
is($msc_2->SO_term, $msc_2_expected, 'vf - most_severe_OverlapConsequence - same rank');

done_testing();
