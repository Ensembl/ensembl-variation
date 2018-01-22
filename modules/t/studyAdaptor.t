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

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');

my $vdb = $multi->get_DBAdaptor('variation');

my $sta = $vdb->get_StudyAdaptor();

ok($sta && $sta->isa('Bio::EnsEMBL::Variation::DBSQL::StudyAdaptor'), "isa study adaptor");


# Values
my $source_name  = 'DGVa';
my $name         = 'estd1';
my %study_list   = ( $name => 4237, 'estd55' => 4246 );
my @study_IDs    = sort {$a <=> $b} values(%study_list);
my $description  = 'Redon 2006 "Global variation in copy number in the human genome." PMID:17122850 [remapped from build NCBI35]';
my $url          = 'ftp://ftp.ebi.ac.uk/pub/databases/dgva/estd1_Redon_et_al_2006';
my $type         = 'Control Set';
my $external_ref = 'pubmed/17122850';


# test fetch by dbID
my $study = $sta->fetch_by_dbID($study_IDs[0]);
ok($study->name() eq $name, "name");
ok($study->description() eq $description, "description");
ok($study->url() eq $url, "url");
ok($study->external_reference() eq $external_ref, "reference");
ok($study->type() eq $type, "type");

# test fetch by name
my $study2 = $sta->fetch_by_name($name);
ok($study2->name() eq $name, "study by name");

# test fetch all by source
my $studies = $sta->fetch_all_by_source($source_name);
ok($studies->[0]->name() eq $name, "study by source");

# test fetch all by dbID list
my $studies2 = $sta->fetch_all_by_dbID_list(\@study_IDs);
ok($studies2->[0]->name() eq $name, "study by dbID list");

# test fetch all by external reference
my $studies3 = $sta->fetch_all_by_external_reference($external_ref);
ok($studies3->[0]->name() eq $name, "study by external reference");

# store

print "\n# Test - store\n";

my $sa = $vdb->get_SourceAdaptor();

my $source = Bio::EnsEMBL::Variation::Source->new( -data_types => ['variation'] );
$source->name( 'test study source' );

ok($sa->store($source), "store");

$source = $sa->fetch_by_name('test study source');

ok($source && $source->name eq 'test study source', "fetch stored source");

my $study3 = Bio::EnsEMBL::Variation::Study->new();
$study3->name('test study');
$study3->source( $source );
$study3->description( 'test study for study->store method' );

ok($sta->store( $study3 ), "store");

my $st = $sta->fetch_by_name('test study');

ok($st && $st->name eq 'test study', "fetch stored study");

done_testing();
