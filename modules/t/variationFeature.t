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
ok($vf->ref_allele_string() eq 'A',     "get ref allele");
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

# test getter/setters

my $v2 = Bio::EnsEMBL::Variation::Variation->new(-name => 'rs12311',
                                                 -source => $source);

ok(test_getter_setter($vf, 'variation', $v2),            'set new variation object');
ok(test_getter_setter($vf, 'map_weight', 4),             'set new map_weight');
ok(test_getter_setter($vf, 'allele_string', 'T/G'),      'set new allele_string');
ok(test_getter_setter($vf, 'variation_name', 'rs21431'), 'set new variation name');
ok(test_getter_setter($vf, 'flank_match', '1'),          'set new flank_match');


# test ambiguity code
ok($vf->ambig_code eq 'W', "ambiguity code");

# test variation class
ok($vf->var_class eq 'SNP', "var class");

# test convert to SNP
#ok($vf->convert_to_SNP, 'convert to SNP'); # Need the ensembl-external repository


# test get all VariationSets
my $var3 = $va->fetch_by_name('rs80359159');
my $vf3 = $var3->get_all_VariationFeatures()->[0];
my $vss = $vf3->get_all_VariationSets();
my $vs_name = 'clinically associated';
ok($vss->[0]->name eq $vs_name, 'get all VariationSets');

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


#test deprecated methods
print "\n## Test deprecated methods ##\n";

# test source
ok($vf->source eq 'dbSNP', "deprecated 'sources'");
# test add consequence type
my $oc2 = Bio::EnsEMBL::Variation::OverlapConsequence->new(-SO_term => 'synonymous_variant');
ok($vf->add_consequence_type($oc2), "deprecated 'add_consequence_type'");
# get consequence type
ok($vf->get_consequence_type()->[0] eq $consequence->SO_term, "deprecated 'get_consequence_type'");


done_testing();
