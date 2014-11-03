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
use Bio::EnsEMBL::Test::MultiTestDB;

use Bio::EnsEMBL::Variation::Variation;
use Bio::EnsEMBL::Variation::Phenotype;
use Bio::EnsEMBL::Variation::PhenotypeFeature;
use Bio::EnsEMBL::Slice;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');

my $vdb = $multi->get_DBAdaptor('variation');
my $db  = $multi->get_DBAdaptor('core');



my $variation = Bio::EnsEMBL::Variation::Variation->new(-name   => 'rs142276873',
                                                        -source => 'dbSNP');



my $sa = $db->get_SliceAdaptor();
my $slice = $sa->fetch_by_region('chromosome', '18');




my $external_id = 12345;
my $p_value     = 0.0000023;
my $risk_allele = 'G';
my $desc        = 'Tea Consumption';
my $gene        = 'TEA1';

my $phenotype = Bio::EnsEMBL::Variation::Phenotype->new(-DESCRIPTION => $desc);

my $pf = Bio::EnsEMBL::Variation::PhenotypeFeature->new(
    -slice     => $slice,
    -start     => 23821095,
    -end       => 23821095,
    -phenotype => $phenotype,
    -type      => 'Variation',
    -object    => $variation,
    -source    => 'OMIM',
    -attribs   => {
      p_value         => $p_value,
      external_id     => $external_id,
      risk_allele     => $risk_allele,
      associated_gene => $gene
    },
    );


ok($pf->start() == 23821095,                   "start");
ok($pf->end()   == 23821095,                   "end") ;
ok($pf->external_id() eq $external_id,         "external_id");
ok($pf->risk_allele() eq $risk_allele,         "risk_allele");
ok($pf->p_value ()    eq  $p_value,            "p_value");
ok($pf->type()        eq 'Variation',          "type");
ok($pf->associated_gene  eq $gene,             "associated_gene");
ok($pf->phenotype()->description eq $desc,     "phenotype");
ok($pf->object()->name()   eq 'rs142276873',   "variation name");


done_testing();
