# Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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
	use Test::More;
	plan tests => 22;
}

use Bio::EnsEMBL::Test::TestUtils;
use_ok('Bio::EnsEMBL::Variation::Variation');
use_ok('Bio::EnsEMBL::Variation::VariationFeature');




# test constructor

my $v = Bio::EnsEMBL::Variation::Variation->new(-name => 'rs2421',
                                                -source => 'dbSNP');


my $start = 100;
my $end   = 100;
my $strand = 1;
my $vname = $v->name();
my $map_weight = 1;
my $allele_str = 'A/T';
my $source = 'dbSNP';
my $is_somatic = 0;
my $minor_allele = 'A';
my $minor_allele_frequency = 0.1;
my $minor_allele_count = 10;

my $vf = Bio::EnsEMBL::Variation::VariationFeature->new
  (-start => 100,
   -end   => 100,
   -strand => 1,
   -variation_name => $vname,
   -map_weight => $map_weight,
   -allele_string => $allele_str,
   -variation => $v,
   -source => $source,
   -is_somatic => $is_somatic,
   -minor_allele  => $minor_allele ,
   -minor_allele_frequency  => $minor_allele_frequency ,
   -minor_allele_count => $minor_allele_count,
);


ok($vf->start()   == $start,            "get start");
ok($vf->end()     == $end,              "get end");
ok($vf->strand()  == $strand,           "get strand");
ok($vf->variation_name() eq $vname,     "get name");
ok($vf->map_weight() == $map_weight,    "get map_weight");
ok($vf->allele_string() eq $allele_str, "get allele");
ok($vf->display_id() eq $vname,         "display_name");
ok($vf->source()    eq $source,         "source");
ok($vf->length()    == 1,               "length");
ok($vf->is_somatic() eq $is_somatic,    "is_somatic");
ok($vf->minor_allele()  eq $minor_allele, "minor allele"); 
ok($vf->minor_allele_frequency() == $minor_allele_frequency , "minor allele freq"); 
ok($vf->minor_allele_count() == $minor_allele_count,  "minor allele count"); 


# test getter/setters

my $v2 = Bio::EnsEMBL::Variation::Variation->new(-name => 'rs12311',
                                                 -source => 'dbSNP');

ok(test_getter_setter($vf, 'variation', $v2));
ok(test_getter_setter($vf, 'map_weight', 4));
ok(test_getter_setter($vf, 'allele_string', 'T|G'));
ok(test_getter_setter($vf, 'variation_name', 'rs21431'));
ok(test_getter_setter($vf, 'flank_match', '1'));


#test ambiguity code
ok($vf->ambig_code eq 'W', "ambiguity code");

#test variation class
ok($vf->var_class eq 'SNP', "class");

