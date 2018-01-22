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

use Bio::EnsEMBL::Variation::Source;
use Bio::EnsEMBL::Variation::Study;
use_ok('Bio::EnsEMBL::Variation::StructuralVariation');
use_ok('Bio::EnsEMBL::Variation::StructuralVariationFeature');


my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');

my $vdb = $multi->get_DBAdaptor('variation');
my $cdb = $multi->get_DBAdaptor('core');

my $svf_adaptor   = $vdb->get_StructuralVariationFeatureAdaptor;
my $slice_adaptor = $cdb->get_SliceAdaptor;

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
ok($svf->bound_start() == $outer_start, "svf -> bound_start");
ok($svf->bound_end()   == $outer_end,   "svf -> bound_end");
ok($svf->variation_name() eq $var_name, "svf -> varname" );
ok($svf->study->name() eq $study->name, "svf -> study" );
ok($svf->is_somatic() eq $is_somatic,   "svf -> is_somatic");
ok($svf->class_SO_term() eq $SO_term,   "svf -> class");
ok($svf->display_id() eq $var_name,     "svf -> display_id");
# source
ok($svf->source->name() eq $source_name,     'svf -> source' );
ok($svf->source_name eq $source_name,               'svf -> source_name');
ok($svf->source_description eq $source_description, 'svf -> source_description');
ok($svf->source_version eq $source_version,         'svf -> source_version');


# test getter/setters
my $var_name2 = 'esv89107';
my $v2 = Bio::EnsEMBL::Variation::StructuralVariation->new(-name   => $var_name2,
                                                           -source => $source);

ok(test_getter_setter($svf, 'structural_variation', $v2), "get/set structural variation");
ok(test_getter_setter($svf, 'variation_name', $var_name2), "get/set name");



is_deeply(
  Bio::EnsEMBL::Variation::StructuralVariationFeature->new_fast({
    class_SO_term => 'deletion',
    start => 11,
    end => 20,
    chr => 1,
  })->to_VCF_record(),
  [1, 10, '.', 'N', '<DEL>', '.', '.', 'END=20'],
  'to_VCF_record'
);


## Other ##

my $svf_id = 4509635;
my $svf2 = $svf_adaptor->fetch_by_dbID($svf_id);

# test get all VariationSets
my $set = '1000 Genomes - High coverage - Trios';
my $vss = $svf2->get_all_VariationSets();
ok($vss->[0]->name eq $set, 'get_all_VariationSets');

# test get nearest Gene
my $gene_stable_id = 'ENSG00000254955';
my $gene = $svf2->get_nearest_Gene();
ok($gene->[0]->stable_id eq $gene_stable_id, 'get_nearest_Gene');

# test get all TranscriptStructuralVariations
my $tr_stable_id = 'ENST00000438775';
my $tsvs = $svf2->get_all_TranscriptStructuralVariations();
ok($tsvs->[0]->transcript->stable_id eq $tr_stable_id, 'get_all_TranscriptStructuralVariations');

# test get all supporting evidence classes
my $secs = $svf2->get_all_supporting_evidence_classes();
ok($secs->[0] eq 'copy_number_loss', 'get_all_supporting_evidence_classes');

# test source object
my $sv_source = $svf2->source();
ok($svf2->source($sv_source), 'source (using argument)');

# test study object
my $sv_study = $svf2->study();
ok($svf2->study($sv_study), 'study object (using argument)');

# test get reference sequence
my $length = $svf2->seq_region_end - $svf2->seq_region_start + 1;
my $ref_seq_length = length($svf2->get_reference_sequence());
ok($ref_seq_length == $length, 'get_reference_sequence');

# test transform
my $svf3 = $svf_adaptor->fetch_by_dbID($svf_id);
my $svf3_chr = $svf3->seq_region_name;
my $svf3_end = $svf3->seq_region_end;
my $contig = 'AC087763.10';
my $svf3_contig = $svf3->transform('contig');
ok($svf3_contig->seq_region_name eq $contig, 'transform to contig');

# test transfer
my $chr_start = 8790818;
my $chr_end   = $svf3_end;
my $slice = $slice_adaptor->fetch_by_region('chromosome', $svf3_chr);
my $svf3_new = $svf3_contig->transfer($slice);
ok($svf3_new->seq_region_name eq $svf3_chr && $svf3_new->seq_region_start == $chr_start && $svf3_new->seq_region_end == $chr_end, 'transfert from contig to chr');


done_testing();

