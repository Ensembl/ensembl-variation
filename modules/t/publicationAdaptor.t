# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2019] EMBL-European Bioinformatics Institute
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

our $verbose = 0;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');

my $vdb = $multi->get_DBAdaptor('variation');

my $va = $vdb->get_VariationAdaptor();
my $pa = $vdb->get_PublicationAdaptor();

ok($pa && $pa->isa('Bio::EnsEMBL::Variation::DBSQL::PublicationAdaptor'), "isa publication adaptor");

# test fetch by PMID

my $pub = $pa->fetch_by_pmid('22779046');


ok($pub->pmid() eq '22779046',       'PMID by PMID');
ok($pub->pmcid() eq 'PMC3392070',   'PMCID by PMID'  );
ok($pub->year() eq '2012',          'year by PMID');
ok($pub->doi() eq '10.1001/2012.journal.123', 'doi by PMID');
ok($pub->ucsc_id() eq 'PMC3392070',   'UCSC ID by PMID'  );
ok($pub->title() eq 'Coanalysis of GWAS with eQTLs reveals disease-tissue associations.',    'title by PMID');
ok($pub->authors() eq 'Kang HP, Morgan AA, Chen R, Schadt EE, Butte AJ.',     'authors by PMID');
ok($pub->variations()->[0]->name() eq 'rs7698608', 'variations by PMID - name of first');
ok(scalar @{$pub->variations()} ==1,            'variations by PMID - count');

my $pub2 = $pa->fetch_by_pmcid('PMC3392070',);
ok($pub2->pmid() eq '22779046',       'PMID by PMCID');

my $pub3 = $pa->fetch_by_dbID(36249);
ok($pub3->pmid() eq '22779046',       'PMID by dbID');

my $pub4 = $pa->fetch_by_doi('10.1001/2012.journal.123');
ok($pub4->pmid() eq '22779046',       'PMID by doi');

my $pubs = $pa->fetch_all_by_dbID_list([36249]);
ok($pubs->[0]->pmid() eq '22779046',       'PMID by dbID list');


my $var = $va->fetch_by_name("rs7698608");
my $pubs2 = $pa->fetch_all_by_Variation($var);
ok($pubs2->[0]->pmid() eq '22779046',     'PMID by variation');
ok($pubs2->[0]->pmcid() eq 'PMC3392070',   'PMCID by variation'  );


## count +/- fails
my $vars = $va->fetch_all_by_publication($pubs->[0]);
ok(scalar(@{$vars}) ==1,   "variation count by publication - no fails");


$va->db->include_failed_variations(1);
my $varfs = $va->fetch_all_by_publication($pubs->[0]);
ok(scalar(@{$varfs}) ==2,   "variation count by publication - inc fails");


## store
my $pub_store = Bio::EnsEMBL::Variation::Publication->new( 
                -title    => "title",
                -authors  => "authorString",
                -pmid     => 1234,
                -pmcid    => "PMC1234",
                -ucsc_id  => "12345",
                -year     => 2010,
                -doi      => "doi:1234",
                -adaptor  => $pa
                );

$pa->store($pub_store);

## try updates
ok($pub_store->ucsc_id() eq '12345',   'inserted UCSC id'  );
$pa->update_ucsc_id($pub_store, "updated");
ok($pub_store->ucsc_id() eq 'updated',   'update UCSC id'  );

## and from db
my $pubup = $pa->fetch_by_dbID($pub_store->dbID);
ok($pubup->ucsc_id() eq 'updated',   'update UCSC id in db');

## update citation
$var = $va->fetch_by_dbID(4770800);
$pa->update_variant_citation($pubup, 615, [$var]);
ok($pubup->variations()->[0]->name() eq 'rs7569578', "citation update");

# Update data source
my $var_2 = $va->fetch_by_dbID(50478027);
$pa->update_variant_citation($pubup, 617, [$var_2]);
$pa->update_citation_data_source(615, 50478027, 1234);

## get publication source
my $pub_1 = $pa->fetch_by_pmid(1234);
$pub_1->get_sources();

my $pub_2 = $pa->fetch_by_pmid(12007216);
$pub_2->get_sources();

my $sources_1 = $pa->fetch_sources_by_pmid(1234);
my $sources_2 = $pa->fetch_sources_by_pmid(12007216);
print Dumper($sources_1), "\n";
print Dumper($sources_2), "\n";

my $sources_3 = $pa->fetch_sources_variation(36250,50478027);
print Dumper($sources_3), "\n";

my $sources_4 = $pa->fetch_sources_variation(20404,26469702);
print Dumper($sources_4), "\n";

done_testing();
