# Copyright [1999-2013] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

BEGIN { $| = 1;
        use Test;
        plan tests => 3;
}


use Bio::EnsEMBL::Test::TestUtils;


use Bio::EnsEMBL::Test::MultiTestDB;

our $verbose = 1;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();

my $vdb = $multi->get_DBAdaptor('variation');
my $db  = $multi->get_DBAdaptor('core');

$vdb->dnadb($db);

my $va  = $vdb->get_VariationAdaptor();
my $iga = $vdb->get_IndividualGenotypeAdaptor();

ok($va && $va->isa('Bio::EnsEMBL::Variation::DBSQL::VariationAdaptor'));
ok($iga && $iga->isa('Bio::EnsEMBL::Variation::DBSQL::CompressedGenotypeAdaptor'));

my $var = $va->fetch_by_dbID(1199);
ok($var->name() eq 'rs1205');

my $igs = $iga->fetch_all_by_Variation($var);

#use Data::Dumper;
#print STDERR Dumper($igs);

print_feats($igs);

sub print_feats {
  my $feats = shift;
  return if(!$verbose);

  foreach my $f (@$feats) {
    print STDERR $f->start(), '-', $f->end(), '-', $f->allele1, '-', $f->allele2, ' ', $f->strand(), "\n";
  }	
	
}	
