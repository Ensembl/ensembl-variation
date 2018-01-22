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



use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Test::MultiTestDB;


our $verbose = 0;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');

my $vdb = $multi->get_DBAdaptor('variation');

my $sa = $vdb->get_SourceAdaptor();

ok($sa && $sa->isa('Bio::EnsEMBL::Variation::DBSQL::SourceAdaptor'), "isa source adaptor");


# Values

my $name           = 'dbSNP';
my $version        = 138;
my $description    = 'Variants (including SNPs and indels) imported from dbSNP';
my $url            = 'http://www.ncbi.nlm.nih.gov/projects/SNP/';
my $type           = undef;
my $somatic_status = 'mixed';
my @data_types     = ('variation');
my %source_list    = ( $name => 1, 'COSMIC' => 26 );
my @source_IDs     = sort {$a <=> $b} values(%source_list);


# test fetch by dbID
my $source = $sa->fetch_by_dbID(1);

ok($source, "source by dbID");
ok($source->name() eq $name, "name");
ok($source->version() == $version, "version");
ok($source->description() eq $description, "$description");
ok($source->url() eq $url, "url");
ok(!$source->type(), "type");
ok($source->somatic_status() eq $somatic_status, "somatic_status");
ok($source->get_all_data_types()->[0] eq $data_types[0], "data_type");

# test fetch by name
my $source2 = $sa->fetch_by_name($name);
ok($source2->name() eq $name, "source by name");

# test fetch all by dbID list
my $sources = $sa->fetch_all_by_dbID_list(\@source_IDs);
ok($source_list{$sources->[0]->name}, 'source by dbID list - 1');
ok($source_list{$sources->[1]->name}, 'source by dbID list - 2');

# test get source version
ok($sa->get_source_version($name) eq $version, "get_source_version");


# store
print "\n# Test - store\n";

delete $source->{$_} for qw(dbID name);
$source->name('test');

ok($sa->store($source), "store");

$source = $sa->fetch_by_name('test');
ok($source && $source->name eq 'test', "fetch stored");

$source->version(140);
$sa->update_version($source);

my $new_source = $sa->fetch_by_name($source->name());
ok($new_source->version() eq 140, "update version");

done_testing();

