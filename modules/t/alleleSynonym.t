# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2023] EMBL-European Bioinformatics Institute
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

use Test::Exception;
use Test::More;

use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Variation::AlleleSynonym;

our $verbose = 0;

# test constructor - with variation object
my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');
my $vdb = $multi->get_DBAdaptor('variation');

my $va  = $vdb->get_VariationAdaptor();
my $asa = $vdb->get_AlleleSynonymAdaptor();

# test constructor
my $as_name = 'CA6124251';
my $hgvs_genomic = 'NC_000011.9:g.66318972G>A';
my $var_name = 'rs490998';
my $var_id = 338639;

my $var = $va->fetch_by_name($var_name);

my $allele_synonym = Bio::EnsEMBL::Variation::AlleleSynonym->new
                (-hgvs_genomic => $hgvs_genomic,
                 -name         => $as_name,
                 -variation    => $var);
                 
ok($allele_synonym->name() eq $as_name, 'as1 - allele synonym name');
ok($allele_synonym->hgvs_genomic() eq  $hgvs_genomic, 'as1 - hgvs_genomic');
ok($allele_synonym->variation()->name() eq $var_name, 'as1 - variation name');

# test getters/test_getter_setters                 
ok(test_getter_setter($allele_synonym, 'name', 'new name'), 'as1 - get/set name');
ok(test_getter_setter($allele_synonym, 'hgvs_genomic', 'new_hgvs'), 'as1 - get/set hgvs');

# test constructor - with $variation_id
my $allele_synonym_2  = Bio::EnsEMBL::Variation::AlleleSynonym->new(
        -_variation_id   => $var_id,
        -adaptor         => $asa,
        -hgvs_genomic    => $hgvs_genomic,
        -name            => $as_name
    );

ok($allele_synonym_2->name() eq $as_name, 'as2 - allele synonym name');
ok($allele_synonym_2->hgvs_genomic() eq  $hgvs_genomic, 'as2 - hgvs_genomic');
ok($allele_synonym_2->variation()->name() eq $var_name, 'as2 - variation name');

# test getters/test_getter_setters                 
ok(test_getter_setter($allele_synonym_2, 'name', 'new name 2'), 'as2 - get/set name');
ok(test_getter_setter($allele_synonym_2, 'hgvs_genomic', 'new_hgvs_2'), 'as2 - get/set hgvs');

# Check 
throws_ok { $allele_synonym_2->variation('test_str'); } 
        qr/Bio::EnsEMBL::Variation::Variation argument expected/, 
        'variation - Throw on invalid variation argument';

done_testing();
