# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2026] EMBL-European Bioinformatics Institute
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
use Test::Deep;
use Test::Exception;

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

## Authors stored as NULL when value not provided ("" or 0)
my $pub_store_null = Bio::EnsEMBL::Variation::Publication->new( 
                -title    => "ABCD",
                -authors  => "",
                -pmid     => 57,
                -pmcid    => "PMC57",
                -ucsc_id  => "12345",
                -year     => 2020,
                -doi      => "doi:12345",
                -adaptor  => $pa
                );

$pa->store($pub_store_null);
my $publication = $pa->fetch_by_dbID($pub_store_null->dbID);
ok(!defined $publication->authors(), "authors NULL");

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

# publication sources
throws_ok { $pub_store->set_variation_id_to_source('var', 'EPMC'); } qr/variation_id is not valid/, 'Throw OK if variation id not valid';
throws_ok { $pub_store->set_variation_id_to_source(4770800, 'Cosmic'); } qr/Source is not valid/, 'Throw OK if publication source not valid';
$pub_store->set_variation_id_to_source(4770800, 'UCSC');
$pub_store->set_variation_id_to_source(4770800, 'dbSNP');

throws_ok { $pub_store->get_all_sources_by_Variation(); } qr/Variation argument is required/, 'Throw OK if variation object is not valid';
my $publication_var_sources = $pub_store->get_all_sources_by_Variation($var);
ok(scalar @$publication_var_sources == 2, 'get_all_sources_by_Variation - Number of sources');

my $expected_source_1 = { 'dbSNP' => 1 };
my $expected_source_2 = { 'EPMC' => 1,
                          'dbSNP' => 1 };

my $variation = $va->fetch_by_dbID(26469702);
my $publications = $pa->fetch_all_by_Variation($variation);
is_deeply($publications->[0]->{variation_id_to_source}->{26469702}, $expected_source_1, "fetch_all_by_Variation - publication sources 1");
is_deeply($publications->[1]->{variation_id_to_source}->{26469702}, $expected_source_2, "fetch_all_by_Variation - publication sources 2");
# Citation without source
my $variation_no_citation = $va->fetch_by_dbID(39404961);
my $publications_no_citation = $pa->fetch_all_by_Variation($variation_no_citation);
ok(!$publications_no_citation->[0]->{variation_id_to_source}, 'fetch_all_by_Variation - source not defined');

# Delete one (and only) publication for a failed variant
my $var_id = 14128071;
$variation = $va->fetch_by_dbID($var_id);
ok($variation->is_failed, "failed variant");
ok(grep(/^Cited$/, @{ $va->fetch_by_dbID($var_id)->get_all_evidence_values() }),
   "variant has 'Cited' evidence");

$publications = $variation->get_all_Publications();
ok(scalar @$publications == 1, "one publication associated with failed variant");
$publications->[0]->variations(); # get all variations

my $vf = $variation->get_all_VariationFeatures();
ok(scalar @$vf == 1, "one variationFeature associated with variant");
ok($vf->[0]->{'display'} == 1, "variationFeature with display = 1");

$pa->remove_publication_by_dbID($publications->[0]->dbID);
ok(!defined $pa->fetch_by_dbID($publications->[0]->dbID), "publication removed");
ok(scalar @{ $variation->get_all_Publications()} == 0,
   "no publications associated with failed variant after removing publication");
ok(!grep(/^Cited$/, @{ $va->fetch_by_dbID($var_id)->get_all_evidence_values() }),
   "No 'Cited' attribute");
ok($variation->get_all_VariationFeatures()->[0]->{'display'} == 0,
   "variationFeature with display = 0 after removing publication");

# Re-add publication and confirm if working
$pa->store($pub, 1);
$pubs = $va->fetch_by_name('rs7698608')->get_all_Publications();
ok($pub->title() eq $pubs->[0]->title(), 'publication re-added with associated variants');

done_testing();
