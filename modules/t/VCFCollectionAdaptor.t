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
use Data::Dumper;
use JSON;


use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Test::MultiTestDB;


our $verbose = 0;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');

my $vdb = $multi->get_DBAdaptor('variation');

# set the VCFCollection config
my $dir = $multi->curr_dir();
ok($vdb->vcf_config_file($dir.'/vcf_config.json') eq $dir.'/vcf_config.json', "DBAdaptor vcf_config_file");
my $vca = $vdb->get_VCFCollectionAdaptor();

ok($vca && $vca->isa('Bio::EnsEMBL::Variation::DBSQL::VCFCollectionAdaptor'), "isa VCFCollectionAdaptor");

# fetch all
my $collections = $vca->fetch_all();
ok($collections && scalar @$collections == 3, "fetch_all count");

# fetch by ID
my $coll = $vca->fetch_by_id('1000genomes_phase1');
ok($coll && $coll->isa('Bio::EnsEMBL::Variation::VCFCollection'), "fetch_by_id isa VCFCollection");

# remove
ok($vca->remove_VCFCollection_by_ID('1000genomes_phase1'), "remove_VCFCollection_by_ID");
is($vca->fetch_by_id('1000genomes_phase1'), undef, 'fetch_by_id after remove');

ok($vca->add_VCFCollection($coll), 'add_VCFCollection');
is($vca->fetch_by_id('1000genomes_phase1') + 0, $coll + 0, 'fetch_by_id after add');

## disallow duplicate IDs
open IN, $dir.'/vcf_config.json';
local $/ = undef;
my $json_string = <IN>;
close IN;
my $config = JSON->new->decode($json_string);


# create duplicate
push @{$config->{collections}}, $config->{collections}->[0];

delete($vca->{collections});

delete $vdb->{vcf_config};
$vdb->vcf_config($config);

throws_ok {$vca->new($vdb)} qr/Collection with ID .+ already exists/, "duplicate collection ID";

done_testing();
