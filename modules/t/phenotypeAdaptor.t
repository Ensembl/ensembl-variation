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
use FindBin qw($Bin);

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::Population;
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Test::MultiTestDB;
use FileHandle;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');
my $vdba = $multi->get_DBAdaptor('variation');

my $pa = $vdba->get_PhenotypeAdaptor();
ok($pa && $pa->isa('Bio::EnsEMBL::Variation::DBSQL::PhenotypeAdaptor'), "get phenotype adaptor");

my $p = $pa->fetch_by_dbID(1);
ok($p && $p->name eq 'ACH', "fetch_by_dbID");

$p = $pa->fetch_by_description('ACHONDROPLASIA')->[0];
ok($p && $p->name eq 'ACH', "fetch_by_description");

$p->name('test');
$p->description('test');
delete $p->{dbID};

ok($pa->store($p), "store");
$p = $pa->fetch_by_description('test')->[0];
ok($p && $p->name eq 'test', "fetch stored");

done_testing();

