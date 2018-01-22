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
use Test::Exception;


our $verbose = 0;


## SETUP
########

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');

my $vdb = $multi->get_DBAdaptor('variation');
my $cdb = $multi->get_DBAdaptor('core');

my $tha = $vdb->get_TranscriptHaplotypeAdaptor();

ok($tha && $tha->isa('Bio::EnsEMBL::Variation::DBSQL::TranscriptHaplotypeAdaptor'), "isa th adaptor");

my $ta = $cdb->get_TranscriptAdaptor;
my $tr = $ta->fetch_by_stable_id('ENST00000502692');

# set the VCFCollection config
my $dir = $multi->curr_dir();
ok($vdb->vcf_config_file($dir.'/th_vcf_config.json') eq $dir.'/th_vcf_config.json', "DBAdaptor vcf_config_file");
my $vca = $vdb->get_VCFCollectionAdaptor();
my $vcf_coll = $vca->fetch_all->[0];
my $temp = $vcf_coll->filename_template();
$temp =~ s/###t\-root###/$dir/;
$vcf_coll->filename_template($temp);
my ($s) = grep {$_->name =~ /NA18499/} @{$vcf_coll->get_all_Samples};

my $c = $tha->get_TranscriptHaplotypeContainer_by_Transcript($tr);
ok($c->isa('Bio::EnsEMBL::Variation::TranscriptHaplotypeContainer'), "get_TranscriptHaplotypeContainer_by_Transcript");

my ($p) = grep {$_->name =~ /CEU/} @{$c->get_all_Populations};

my $h = $tha->fetch_all_by_Transcript($tr);
is(scalar @$h, 116, "fetch_all_by_Transcript - count");
ok($h->[0]->isa('Bio::EnsEMBL::Variation::TranscriptHaplotype'), "fetch_all_by_Transcript - isa");

$h = $tha->fetch_all_by_Transcript($tr, $s);
is(scalar @$h, 4, "fetch_all_by_Transcript with sample - count");
ok($h->[0]->isa('Bio::EnsEMBL::Variation::TranscriptHaplotype'), "fetch_all_by_Transcript with sample - isa");

$h = $tha->fetch_all_by_Transcript($tr, $p);
is(scalar @$h, 13, "fetch_all_by_Transcript with population - count");
ok($h->[0]->isa('Bio::EnsEMBL::Variation::TranscriptHaplotype'), "fetch_all_by_Transcript with population - isa");

throws_ok(
  sub {$tha->fetch_all_by_Transcript($tr, $c)},
  qr/Expected.+Sample.+Population/,
  "fetch_all_by_Transcript with wrong object type 1"
);


$h = $tha->fetch_all_CDSHaplotypes_by_Transcript($tr);
is(scalar @$h, 75, "fetch_all_CDSHaplotypes_by_Transcript - count");
ok($h->[0]->isa('Bio::EnsEMBL::Variation::CDSHaplotype'), "fetch_all_CDSHaplotypes_by_Transcript - isa");

$h = $tha->fetch_all_CDSHaplotypes_by_Transcript($tr, $s);
is(scalar @$h, 2, "fetch_all_CDSHaplotypes_by_Transcript with sample - count");
ok($h->[0]->isa('Bio::EnsEMBL::Variation::CDSHaplotype'), "fetch_all_CDSHaplotypes_by_Transcript with sample - isa");

$h = $tha->fetch_all_CDSHaplotypes_by_Transcript($tr, $p);
is(scalar @$h, 8, "fetch_all_CDSHaplotypes_by_Transcript with population - count");
ok($h->[0]->isa('Bio::EnsEMBL::Variation::CDSHaplotype'), "fetch_all_CDSHaplotypes_by_Transcript with population - isa");


$h = $tha->fetch_all_ProteinHaplotypes_by_Transcript($tr);
is(scalar @$h, 41, "fetch_all_ProteinHaplotypes_by_Transcript - count");
ok($h->[0]->isa('Bio::EnsEMBL::Variation::ProteinHaplotype'), "fetch_all_ProteinHaplotypes_by_Transcript - isa");

$h = $tha->fetch_all_ProteinHaplotypes_by_Transcript($tr, $s);
is(scalar @$h, 2, "fetch_all_ProteinHaplotypes_by_Transcript with sample - count");
ok($h->[0]->isa('Bio::EnsEMBL::Variation::ProteinHaplotype'), "fetch_all_ProteinHaplotypes_by_Transcript with sample - isa");

$h = $tha->fetch_all_ProteinHaplotypes_by_Transcript($tr, $p);
is(scalar @$h, 5, "fetch_all_ProteinHaplotypes_by_Transcript with population - count");
ok($h->[0]->isa('Bio::EnsEMBL::Variation::ProteinHaplotype'), "fetch_all_ProteinHaplotypes_by_Transcript with population - isa");


# check filtering works
is(scalar @{$c->get_all_TranscriptHaplotypes}, 116, 'get_TranscriptHaplotypeContainer_by_Transcript - no frequency filter');
is(scalar (grep {$_->name =~ /290P\>L/} @{$c->get_all_ProteinHaplotypes}), 1, 'no frequency filter check included');

$c = $tha->get_TranscriptHaplotypeContainer_by_Transcript($tr, {frequency => {}});
is(scalar @{$c->get_all_TranscriptHaplotypes}, 27, 'get_TranscriptHaplotypeContainer_by_Transcript - frequency filter count');
is(scalar (grep {$_->name =~ /290P\>L/} @{$c->get_all_ProteinHaplotypes}), 0, 'frequency filter check excluded 1');
is(scalar (grep {$_->name =~ /254R\>Q/} @{$c->get_all_ProteinHaplotypes}), 0, 'frequency filter check excluded 2');

$c = $tha->get_TranscriptHaplotypeContainer_by_Transcript($tr, {frequency => {frequency => 0.005}});
is(scalar @{$c->get_all_TranscriptHaplotypes}, 37, 'get_TranscriptHaplotypeContainer_by_Transcript - frequency filter count 0.005');
is(scalar (grep {$_->name =~ /254R\>Q/} @{$c->get_all_ProteinHaplotypes}), 1, 'frequency filter 0.05 check included');

$c = $tha->get_TranscriptHaplotypeContainer_by_Transcript($tr, {frequency => {population => '1000GENOMES:phase_3:AFR'}});
is(scalar @{$c->get_all_TranscriptHaplotypes}, 82, 'get_TranscriptHaplotypeContainer_by_Transcript - frequency filter count AFR');

done_testing();
