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
use Bio::EnsEMBL::Test::MultiTestDB;

BEGIN {
    use_ok('Bio::EnsEMBL::Variation::Utils::SeqRegionUtils', qw(update_seq_region_ids));
}

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');
my $cdb = $multi->get_DBAdaptor('core');
my $vdb = $multi->get_DBAdaptor('variation');

## check no seq_region update performed when last mapping_set entry is not related to current release
# NOTE: test databases do not have the release and assembly on dbname
my $dbname = $vdb->dbc->dbname;
my $dbh = $vdb->dbc->db_handle;
my $sthCheck = $dbh->prepare("SELECT UPDATE_TIME
                              FROM   information_schema.tables
                              WHERE  TABLE_SCHEMA = '$dbname'
                              AND TABLE_NAME = 'seq_region';");
$sthCheck->execute();
my ($update_time_before) = $sthCheck->fetchrow_array();
update_seq_region_ids($cdb,$vdb);
$sthCheck->execute();
my ($update_time_after) = $sthCheck->fetchrow_array();
ok($update_time_before eq $update_time_after, "no update");

done_testing();
