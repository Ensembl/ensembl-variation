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
use Data::Dumper;
use Test::Exception;
use Test::More;
use Bio::EnsEMBL::Test::MultiTestDB;

use FindBin qw($Bin);

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Variation::Population;
use Bio::EnsEMBL::Variation::StructuralVariationPopulationFrequency;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');

my $vdb  = $multi->get_DBAdaptor('variation');

my $pop_adaptor  = $vdb->get_PopulationAdaptor;
my $svpf_adaptor = $vdb->get_StructuralVariationPopulationFrequencyAdaptor;


# test constructor

## need Population object 
my $pop_id   = 373522;
my $pop_name = '1000GENOMES:phase_3:GBR';
my $pop_desc = 'British in England and Scotland';
my $pop_size = 91;

my $pop = Bio::EnsEMBL::Variation::Population->new(
   -adaptor     => $pop_adaptor,
   -dbID        => $pop_id,
   -name        => $pop_name,
   -description => $pop_desc,
   -size        => $pop_size
);

my $sv_dbID = 90344156;

my $svpf = Bio::EnsEMBL::Variation::StructuralVariationPopulationFrequency->new(
  -adaptor                  => $svpf_adaptor,
  -_structural_variation_id => $sv_dbID,
  -_population_id           => $pop_id,
  -name                     => $pop_name,
  -description              => $pop_desc,
  -size                     => $pop_size,
  -region_name              => 'X',
  -samples_class            => { 'copy_number_gain' => { '18085' => 2, # 1000GENOMES:phase_3:HG00253
                                                         '18003' => 1  # 1000GENOMES:phase_3:HG00106
                                                       },
                                 'copy_number_loss' => { '18832' => 2, # 1000GENOMES:phase_3:HG02215
                                                         '18609' => 1, # 1000GENOMES:phase_3:HG01789
                                                         '18089' => 1  # 1000GENOMES:phase_3:HG00257
                                                       }
                               }
);

$svpf->population($pop); # Bio::EnsEMBL::Population object

ok($svpf->population->dbID() eq $pop_id, 'dbID');
ok($svpf->name() eq $pop_name, 'population name');
ok($svpf->description() eq $pop_desc, 'population description');
ok($svpf->size() == $pop_size, 'population size');

# Frequency
my $gfreqs = sprintf("%.5f",$svpf->frequency)."\n";
ok($gfreqs == 0.05147, 'count global frequency for this population');

my $freqs = $svpf->frequencies_by_class_SO_term();
is_deeply(
  $freqs, 
  {'copy_number_loss' => 0.0294117647058823529, 'copy_number_gain' => 0.0220588235294117647},
  'compare SO term frequencies'
);

# Allele count
ok($svpf->allele_count == 7 , 'count global allele count for this population');

my $counts = $svpf->allele_count_by_class_SO_term();
ok($counts->{'copy_number_gain'} == 3 , 'count "copy_number_gain" allele count for this population');
ok($counts->{'copy_number_loss'} == 4 , 'count "copy_number_loss" allele count for this population');

done_testing();
