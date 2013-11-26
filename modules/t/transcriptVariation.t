# Copyright [1999-2013] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

BEGIN { $| = 1;
	use Test;
	plan tests => 15;
}


use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Variation::Variation;
use Bio::EnsEMBL::Variation::VariationFeature;
use Bio::EnsEMBL::Variation::TranscriptVariation;
use Bio::EnsEMBL::Transcript;

our $verbose = 5;

# test constructor

my $tr = Bio::EnsEMBL::Transcript->new();

my $v = Bio::EnsEMBL::Variation::Variation->new(-name => 'rs665',
                                                -source => 'dbSNP');

my $vf = Bio::EnsEMBL::Variation::VariationFeature->new
  (-start => 23650516,
   -end   => 23650516,
   -strand => -1,
   -variation_name => 'rs665',
   -map_weight => 1,
   -allele_string => 'G/A',
   -variation => $v,
   -dbID => 582);



my $pep_allele = 'V/E';

my $cdna_start = 775;
my $cdna_end   = 775;
my $tl_start   = 255;
my $tl_end     = 255;
my $consequence_type       = 'NON_SYNONYMOUS_CODING';


my $trvar = Bio::EnsEMBL::Variation::TranscriptVariation->new
  (-variation_feature => $vf,
   -transcript        => $tr,
   -pep_allele_string => $pep_allele,
   -cdna_start        => $cdna_start,
   -cdna_end          => $cdna_end,
   -translation_start => $tl_start,
   -translation_end   => $tl_end,
   -consequence_type  => $consequence_type);

ok($trvar->variation_feature() == $vf);
ok($trvar->transcript() == $tr);
ok($trvar->pep_allele_string() eq $pep_allele);
ok($trvar->cdna_start() == $cdna_start);
ok($trvar->cdna_end() == $cdna_end);
ok($trvar->translation_start() == $tl_start);
ok($trvar->translation_end() == $tl_end);
ok($trvar->consequence_type() eq $consequence_type);


# test getter/setters
my $tr_new = Bio::EnsEMBL::Transcript->new();
ok(test_getter_setter($trvar, 'transcript', $tr_new));




ok(test_getter_setter($trvar, 'pep_allele_string', $pep_allele));

ok(test_getter_setter($trvar, 'cdna_start', 1));
ok(test_getter_setter($trvar, 'cdna_end', 12));

ok(test_getter_setter($trvar, 'translation_start', 4));
ok(test_getter_setter($trvar, 'translation_end', 10));

ok(test_getter_setter($trvar, 'consequence_type', 'INTRONIC'));


