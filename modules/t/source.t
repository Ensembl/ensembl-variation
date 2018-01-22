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
use Bio::EnsEMBL::Variation::Source;


my $name           = 'dbSNP';
my $version        = 138;
my $description    = 'Variants (including SNPs and indels) imported from dbSNP (mapped to GRCh38)';
my $url            = 'http://www.ncbi.nlm.nih.gov/projects/SNP/';
my $type           = 'test';
my $somatic_status = 'mixed';
my @data_types     = ('variation','variation_synonym');



# test constructor
my $source = Bio::EnsEMBL::Variation::Source->new
  (-name           => $name,
   -version        => $version,
   -description    => $description,
   -url            => $url,
   -type           => $type,
   -somatic_status => $somatic_status,
   -data_types     => \@data_types
);


ok($source->name() eq $name, "name");
ok($source->version() == $version, "version");
ok($source->description() eq $description, "$description");
ok($source->url() eq $url, "url");
ok($source->type() eq $type, "type");
ok($source->somatic_status() eq $somatic_status, "somatic_status");
ok($source->get_all_data_types()->[0] eq $data_types[0], "data_types 1");
ok($source->get_all_data_types()->[1] eq $data_types[1], "data_types 2");


# test getter/setters
ok(test_getter_setter($source, 'name', 'new name'), "get/set name");
ok(test_getter_setter($source, 'version', 141), "get/set version");
ok(test_getter_setter($source, 'description', 'new description'), "get/set description");
ok(test_getter_setter($source, 'url', 'http://www.ensembl.org'), "get/set url");

my %versions = ($version => $version, '20152' => '2015.2', '201502' => '02/2015', '20150212' => '12/02/2015');
my $count = 1;
foreach my $v (keys(%versions)) {
  $source->version($v);
  ok($source->formatted_version() eq $versions{$v}, "formatted_version - $count");
  $count++;
}


done_testing();
