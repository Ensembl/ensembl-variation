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
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Test::MultiTestDB;

use Bio::EnsEMBL::Variation::Source;
use Bio::EnsEMBL::Variation::Study;
use_ok('Bio::EnsEMBL::Variation::StructuralVariation');
use_ok('Bio::EnsEMBL::Variation::StructuralVariationFeature');


my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');

my $vdb  = $multi->get_DBAdaptor('variation');

my $svf_adaptor = $vdb->get_StructuralVariationFeatureAdaptor;


# test constructor

## need Source object 
my $source_name           = 'DGVa';
my $source_version        = 201310;
my $source_description    = 'Database of Genomic Variants Archive';
my $source_id             = 11;

my $source = Bio::EnsEMBL::Variation::Source->new
  (-dbID        => $source_id,
   -name        => $source_name,
   -version     => $source_version,
   -description => $source_description
);

## need Study object 
my $study_name = 'estd59';
my $study_url = 'ftp://ftp.ebi.ac.uk/pub/databases/dgva/estd59_1000_Genomes_Consortium_Pilot_Project';
my $study_description = '1000 Genomes Project Consortium - Pilot Project. PMID:20981092';

my $study = Bio::EnsEMBL::Variation::Study->new
  (-name         => $study_name,
   -url          => $study_url,
   -description  => $study_description,
   -_source_id   => $source->dbID
);


my $dbID = 4509635;
my $outer_start = 7803891;
my $inner_start = 7805991;
my $inner_end = 7823440;
my $outer_end = 7825340;
my $var_name = 'esv93078';
my $chr = '8';
my $is_somatic = 0;
my $sv_length = $outer_end-$outer_start+1;
my $SO_term = 'copy_number_variant';

## need a StructuralVariation object 
my $sv = Bio::EnsEMBL::Variation::StructuralVariation->new(-name   => $var_name,
                                                           -source => $source);

my $svf = Bio::EnsEMBL::Variation::StructuralVariationFeature->new
  (-adaptor     => $svf_adaptor,
   -outer_start => $outer_start,
   -start       => $outer_start,
   -inner_start => $inner_start,
   -inner_end   => $inner_end,
   -end         => $outer_end,
   -outer_end   => $outer_end,
   -strand      => 1,
   -variation_name => $var_name,
   -structural_variation => $sv,
   -source => $source,
   -study  => $study,
   -is_somatic => $is_somatic,
   -length => $sv_length,
   -class_SO_term => $SO_term,
);


ok($svf->outer_start() == $outer_start, "svf -> outer_start");
ok($svf->start()       == $outer_start, "svf -> start");
ok($svf->inner_start() == $inner_start, "svf -> inner_start");
ok($svf->inner_end()   == $inner_end,   "svf -> inner_end");
ok($svf->end()         == $outer_end,   "svf -> end") ;
ok($svf->outer_end()   == $outer_end,   "svf -> outer_end");
ok($svf->strand()      == 1,            "svf -> strand");
ok($svf->variation_name() eq $var_name, "svf -> varname" );
ok($svf->source_object->name() eq $source->name, "svf -> source" );
ok($svf->study->name() eq $study->name, "svf -> study" );
ok($svf->is_somatic() eq $is_somatic,   "svf -> is_somatic");
ok($svf->class_SO_term() eq $SO_term,   "svf -> class");


# test getter/setters

my $var_name2 = 'esv89107';
my $v2 = Bio::EnsEMBL::Variation::StructuralVariation->new(-name   => $var_name2,
                                                           -source => $source);

ok(test_getter_setter($svf, 'structural_variation', $v2));
ok(test_getter_setter($svf, 'variation_name', $var_name2));

done_testing();

