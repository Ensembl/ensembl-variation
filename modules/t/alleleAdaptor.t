# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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
use FileHandle;

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

# store method
my $pa = $vdba->get_PopulationAdaptor();

my $al_hash = {
  allele     => 'A',
  frequency  => 0.12345,
  population => $pa->fetch_by_name('1000GENOMES:phase_1_GBR'),
  subsnp     => 12345,
  count      => 12345,
  variation  => $variation
};

my $al1 = Bio::EnsEMBL::Variation::Allele->new_fast($al_hash);

ok($aa->store($al1), "store");

$al_hash->{population} = $pa->fetch_by_name('1000GENOMES:phase_1_LWK');
my $al2 = Bio::EnsEMBL::Variation::Allele->new_fast($al_hash);

my $hash2;
%$hash2 = %$al_hash;
$hash2->{frequency} = 1 - $hash2->{frequency};
my $al3 = Bio::EnsEMBL::Variation::Allele->new_fast($hash2);

ok($aa->store_multiple([$al2, $al3]), "store_multiple");

my $fetched = $variation->get_all_Alleles();

ok((grep {$_->frequency && $_->frequency == 0.12345} @$fetched), "fetch stored allele 1");
ok((grep {$_->frequency && $_->frequency == 0.12345} grep {$_->population && $_->population->name eq '1000GENOMES:phase_1_LWK'} @$fetched), "fetch stored allele 2");
ok((grep {$_->frequency && $_->frequency == 0.87655} grep {$_->population && $_->population->name eq '1000GENOMES:phase_1_LWK'} @$fetched), "fetch stored allele 3");

my $fh = FileHandle->new;
my $tmpfile = "$Bin\/$$\_allele.txt";
$fh->open(">$tmpfile") or die("ERROR: Could not write to $tmpfile\n");
ok($aa->store_to_file_handle($al1, $fh), "store_to_file_handle");
$fh->close();

open FH, $tmpfile;
my @lines = <FH>;
close FH;

unlink($tmpfile);

ok(scalar @lines == 1 && $lines[0] =~ /12345/, "check dumped data");

done_testing();

