# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2017] EMBL-European Bioinformatics Institute
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

my $svpf_adaptor = $vdb->get_StructuralVariationSampleAdaptor;


# test constructor

## need Population object 
my $pop_id   = 373522;
my $pop_name = '1000GENOMES:phase_3:GBR';
my $pop_desc = 'British in England and Scotland';
my $pop_size = 91;

my $pop = Bio::EnsEMBL::Variation::Population->new
  (-dbID        => $pop_id,
   -name        => $pop_name,
   -description => $pop_desc,
   -size        => $pop_size
);

my $dbID = 6107305;

my $svpf = Bio::EnsEMBL::Variation::StructuralVariationPopulationFrequency->new(
  -_population_id     => $pop_id,
  -name               => $pop_name,
  -description        => $pop_desc,
  -size               => $pop_size,
  -samples_class      => { 'copy_number_gain' => { '1000GENOMES:phase_3:HG00253' => 'homozygous', 
                                                   '1000GENOMES:phase_3:HG00106' => 'heterozygous'
                                                 },
                           'copy_number_loss' => { '1000GENOMES:phase_3:HG02215' => 'homozygous', 
                                                   '1000GENOMES:phase_3:HG01789' => 'heterozygous',
                                                   '1000GENOMES:phase_3:HG00257' => 'heterozygous'
                                                 }
                         }
);

$svpf->population($pop); # Bio::EnsEMBL::Population object

ok($svpf->population->dbID() eq $pop_id, 'dbID');
ok($svpf->name() eq $pop_name, 'population name');
ok($svpf->description() eq $pop_desc, 'population description');
ok($svpf->size() == $pop_size, 'population size');
my $gfreqs = sprintf("%.5f",$svpf->freqs)."\n";
ok($gfreqs == 0.03846, 'count global frequency for this population');

my $freqs = $svpf->freqs_by_class_SO_term();
ok(scalar(keys(%$freqs)) == 2, 'count number of SO terms');
my $SO_freq = sprintf("%.5f",$freqs->{'copy_number_loss'});
ok(sprintf("%.5f",$freqs->{'copy_number_loss'}) == 0.02198, 'count frequency for SO term "copy_number_loss"');

done_testing();
