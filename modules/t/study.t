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
use Data::Dumper;

use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Test::MultiTestDB;

use Bio::EnsEMBL::Variation::Study;
use Bio::EnsEMBL::Variation::Source;


my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');

my $vdb = $multi->get_DBAdaptor('variation');


## need source object 
my $source_name = 'DGVa';
my $srca   = $vdb->get_SourceAdaptor();
my $source = $srca->fetch_by_name($source_name);


my $name                = 'studyname';
my $description         = 'longer study description';
my $url                 = 'http://www.ensembl.org/';
my $external_reference  = 'pubmed/22566624';
my $type                = 'GWAS';
my %associated_studies  = ( 'estd1' => 4237, 'estd55' => 4246 ); # Examples from the test database
my @asso_study_IDs      = values(%associated_studies);

# test constructor
my $study = Bio::EnsEMBL::Variation::Study->new(
   -adaptor            => $source->adaptor,
   -name               => $name,
   -description        => $description,
   -url                => $url,
   -external_reference => $external_reference, 
   -type               => $type,
   -_source_id         => $source->dbID(),
   -associate          => \@asso_study_IDs
);

#print Dumper $study;
ok($study->name() eq $name, "name");
ok($study->description() eq $description, "$description");
ok($study->url() eq $url, "url");
ok($study->external_reference() eq $external_reference, "reference");
ok($study->type() eq $type, "type");

my $asso_studies = $study->associated_studies();
ok($associated_studies{$asso_studies->[0]->name}, 'associated_studies - 1');
ok($associated_studies{$asso_studies->[1]->name}, 'associated_studies - 2');

## Source methods ##
print "\n## Source methods ##\n";
ok($study->source()->name eq $source_name, "source");
ok($study->source_name() eq $source_name, "source_name");
ok($study->source_version() eq 201310, "source_version");

# test getter/setters
print "\n## getter/setters ##\n";
ok(test_getter_setter($study, 'name', 'new name'), "get/set name");
ok(test_getter_setter($study, 'description', 'new description'), "get/set description");
ok(test_getter_setter($study, 'url', 'http://www.ebi.ac.uk/ega'), "get/set url");


done_testing();
