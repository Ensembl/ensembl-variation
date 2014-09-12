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
use FindBin qw($Bin);

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::Population;
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Test::MultiTestDB;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');
my $vdba = $multi->get_DBAdaptor('variation');


my $name = 'rs144235347';

my $va = $vdba->get_VariationAdaptor();
my $variation = $va->fetch_by_name($name);
my $aa = $vdba->get_AlleleAdaptor();

my $alleles = $aa->fetch_all_by_Variation($variation);

is(scalar @$alleles, 24, 'alleles for variation');
my $hash;
foreach (@$alleles) {
    my $allele = $_->allele();
    my $population = $_->population();
    my $freq = $_->frequency();
    $freq ||= 'No frequency';
    if (defined $population) {
        $hash->{$population->name}->{$allele} = $freq;
        #print $population->name, ' ', $allele, ' ', $freq, "\n";
    }
}
is($hash->{'1000GENOMES:phase_1_LWK'}->{'T'}, 0.0103092783505155, 'allele freq for population');
is($hash->{'1000GENOMES:phase_1_GBR'}->{'C'}, 1, 'allele freq for population');

done_testing();

