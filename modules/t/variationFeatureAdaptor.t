# Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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
my $db  = $multi->get_DBAdaptor('core');

$vdb->dnadb($db);

my $vfa = $vdb->get_VariationFeatureAdaptor();
my $va  = $vdb->get_VariationAdaptor();

ok($vfa && $vfa->isa('Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor'));

my $sa = $db->get_SliceAdaptor();

my $slice = $sa->fetch_by_region('chromosome', '18');

my $vfs = $vfa->fetch_all_by_Slice($slice);

#print Dumper $vfs;
#my $n = @$vfs;
#print "$n\n";
#ok(@$vfs == 68 , "variationfeature count") ;

my $vf = $vfa->fetch_by_dbID(33303674);

ok($vf->start() == 23821095,                 "vf_id -> start");
ok($vf->end()   == 23821095,                 "vf_id -> end") ;
ok($vf->strand() == 1,                       "vf_id -> strand");
ok($vf->allele_string() eq 'G/A',            "vf_id -> allele_string");
ok($vf->variation()->name() eq 'rs142276873',"vf_id -> varname" );
ok($vf->display_id() eq 'rs142276873',       "vf_id -> display id");
ok($vf->map_weight() == 1,                   "vf_id -> map weight");
ok($vf->slice()->name() eq $slice->name(),   "vf_id -> slice name");


# test fetch_all_by_Variation

my $v = $va->fetch_by_dbID(30220007);
$vfs = $vfa->fetch_all_by_Variation($v);
ok(@$vfs == 1,                              "var -> vf count ");

$vf = $vfs->[0];

ok($vf->dbID() == 33303674,                  "var -> vf id");
ok($vf->slice->name() eq $slice->name(),     "var -> slice name ");
ok($vf->start == 23821095,                   "var -> start");
ok($vf->end() == 23821095,                   "var -> end");
ok($vf->strand() == 1,                       "var -> strand");
ok($vf->allele_string() eq 'G/A',            "var -> allele string");
ok($vf->variation()->name() eq 'rs142276873',"var -> name");
ok($vf->display_id() eq 'rs142276873',       "var -> display id");
ok($vf->map_weight() == 1,                   "var -> map weight");
ok($vf->slice()->name() eq $slice->name(),   "var -> slice name");
ok(@{$vf->consequence_type()} == 2 , "var -> consequence type count"); 
ok($vf->display_consequence() eq 'intron_variant', "var -> display consequence"); 
ok($vf->consequence_type()->[0] eq 'NMD_transcript_variant', "var -> consequence type"); 


my $cons = $vf->most_severe_OverlapConsequence();
ok($cons->description() eq 'A transcript variant occurring within an intron', "consequence description");
ok($cons->label() eq 'Intron variant',                                        "consequence label");
ok($cons->SO_term() eq 'intron_variant',                                      "consequence SO_term"); 
ok($cons->SO_accession() eq 'SO:0001627',                                     "consequence SO_accession"); 
ok($cons->tier() eq '3',                                                      "consequence tier"); 
ok($cons->rank() eq '21',                                                     "consequence rank"); 
ok($cons->NCBI_term() eq 'intron',                                            "consequence NCBI term"); 
ok($cons->impact() eq 'MODIFIER',                                             "consequence impact"); 
ok($cons->display_term() eq 'INTRONIC',                                       "consequence display_term"); 



done_testing();

