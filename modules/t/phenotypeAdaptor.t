# Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

my $omulti = Bio::EnsEMBL::Test::MultiTestDB->new('multi');
my $odb = $omulti->get_DBAdaptor('ontology');
Bio::EnsEMBL::Registry->add_db($omulti, 'ontology', $odb);


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


## check ontology accession handling
$p->add_ontology_accession('Orphanet:15', 'Manual');

ok($p->ontology_accessions()->[0] eq 'Orphanet:15', "get ontology accession");
ok($p->ontology_accessions_with_source()->[0]->{source}    eq 'Manual', "get ontology source");
ok($p->ontology_accessions_with_source()->[0]->{accession} eq 'Orphanet:15', "get ontology accession 2");
$pa->store_ontology_accessions($p);
my $pbo1 = $pa->fetch_all_by_ontology_accession('Orphanet:15');
ok($pbo1->[0]->description() eq 'test', "fetch by newly stored ontology accession");

my $pbo = $pa->fetch_all_by_ontology_accession('Orphanet:130');

ok($pbo->[0]->description() eq 'BRUGADA SYNDROME', "fetch by ontology accession"); 
ok($pbo->[0]->ontology_accessions()->[0] eq 'Orphanet:130', "get ontology accession from db");

ok($pbo->[0]->ontology_accessions_with_source()->[0]->{source}    eq 'OLS exact', "get ontology source from db"); 

## test look up by Ontology Synonym
my $ontterm_ad  = $odb->get_OntologyTermAdaptor;
my $terms = $ontterm_ad->fetch_all_by_name("%Dream disease%", "EFO");

my $p_by_OT = $pa->fetch_by_OntologyTerm( $terms->[0] )->[0];
ok($p_by_OT->description() eq 'BRUGADA SYNDROME', "fetch by OntologyTerm synonym");


done_testing();

