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

use Test::Exception;
use Test::More;

use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Test::MultiTestDB;

our $verbose = 0;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');

my $vdb = $multi->get_DBAdaptor('variation');

my $va = $vdb->get_VariationAdaptor();
my $asa = $vdb->get_AlleleSynonymAdaptor();

ok($asa && $asa->isa('Bio::EnsEMBL::Variation::DBSQL::AlleleSynonymAdaptor'), 
                     'isa AlleleSynonymAdaptor');

# Test fetch_by_dbID
throws_ok { $asa->fetch_by_dbID(); } 
           qr/dbID argument expected/, 
           'fetch_by_dbID - Throw on missing dbID argument';

#   Values
my $as_name = 'CA6124251';
my $hgvs_genomic = 'NC_000011.9:g.66318972G>A';
my $var_name = 'rs490998';

my $allele_synonym = $asa->fetch_by_dbID(120201942);
ok($allele_synonym, 'fetch_by_dbID - allele_synonym');
ok($allele_synonym->name() eq $as_name, 'fetch_by_dbID - name');
ok($allele_synonym->hgvs_genomic() eq $hgvs_genomic, 
          'fetch_by dbID - hgvs_genomic');
ok($allele_synonym->variation()->name() eq $var_name, 
          'fetch_by_dbID - variation');

# Test fetch_all_by_Variation
throws_ok { $asa->fetch_all_by_Variation(); } 
        qr/variation argument is required/, 
        'fetch_all_by_Variation - Throw on missing variation argument';

my $var = $va->fetch_by_name('rs509556');
my $allele_synonyms_2 = $asa->fetch_all_by_Variation($var);
ok(@{$allele_synonyms_2} == 2, 'fetch_all_by_Variation - number');

my @expected_as_names = ('CA13462540','CA475368754');
my @as_names = map {$_->name()} @{$allele_synonyms_2};
is_deeply(\@as_names, \@expected_as_names, 'fetch_all_by_Variation - name');

# Test fetch_all_by_name
throws_ok { $asa->fetch_all_by_name(); } 
        qr/name argument expected/, 
        'fetch_all_by_name - Throw on missing name argument';

my $as_name_2 = 'CA224064043';
my @expected_var_names = ('rs370045702','rs61393902');
my $allele_synonyms_3 = $asa->fetch_all_by_name($as_name_2);
ok(@{$allele_synonyms_3} == 2, 'fetch_all_by_name - number');

my @var_names = map {$_->variation()->name()} @{$allele_synonyms_3};
is_deeply(\@var_names, \@expected_var_names, 'fetch_all_by_name - var name');


done_testing();
