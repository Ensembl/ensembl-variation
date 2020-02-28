# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2020] EMBL-European Bioinformatics Institute
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
use Test::Exception;
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
#        print $population->name, ' ', $allele, ' ', $freq, "\n";

    }
}
is(sprintf("%.6f", $hash->{'1000GENOMES:phase_1_LWK'}->{'T'}), 0.010309, 'allele freq for population');
is($hash->{'1000GENOMES:phase_1_GBR'}->{'C'}, 1, 'allele freq for population');

my $pa = $vdba->get_PopulationAdaptor();
my $population_name = '1000GENOMES:phase_1_AFR';
my $population = $pa->fetch_by_name($population_name);

$alleles = $aa->fetch_all_by_Variation($variation, $population);
is(scalar @$alleles, 2, 'alleles for variation and population');

$aa->{_cache} = undef;
$alleles = $aa->fetch_all_by_Variation($variation, $population);
is(scalar @$alleles, 2, 'alleles for variation and population');

my $allele = $aa->fetch_by_dbID(104844866);
$population = $pa->fetch_by_dbID(650);
my $handle = $aa->get_subsnp_handle($allele, $population);
ok($handle eq 'PERLEGEN', 'get_subsnp_handle');

$allele = $aa->fetch_by_dbID(104845982);
my $descriptions = $aa->get_all_failed_descriptions($allele);
ok($descriptions->[0] eq 'Additional submitted allele data from dbSNP does not agree with the dbSNP refSNP alleles', 'get_all_failed_descriptions');

my $allele_code = $aa->_allele_code('A');
ok($allele_code == 2, '_allele_code');
$aa->db->{_allele_codes} = undef;
$allele_code = $aa->_allele_code('A');
ok($allele_code == 2, '_allele_code');

# store method
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

# fetch_all
throws_ok { $aa->fetch_all; } qr/fetch_all cannot be used for Allele objects/, 'Cannot use fetch_all on AlleleAdaptor';

# fetch_all_by_subsnp_id
my $al4 = $aa->fetch_all_by_subsnp_id('ss24191327');
ok(scalar @$al4 == 24, 'fetch_all_by_subsnp_id');
throws_ok { $aa->fetch_all_by_subsnp_id; } qr/name argument expected/, 'Throw on missing name argument';


# read ESP allele frequencies from VCF
my $dir = $multi->curr_dir();
ok($vdba->vcf_config_file($dir.'/vcf_config.json') eq $dir.'/vcf_config.json', "DBAdaptor vcf_config_file");
my $vca = $vdba->get_VCFCollectionAdaptor();
my $coll = $vca->fetch_by_id('esp_GRCh37');
my $temp = $coll->filename_template();
$temp =~ s/###t\-root###/$dir/;
$coll->filename_template($temp);
$aa->db->use_vcf(1);

$variation = $va->fetch_by_name('rs35099512');
my $al5 = $aa->fetch_all_by_Variation($variation);

is_deeply(
  [
    map {'p:'.$_->population->name.' a:'.$_->allele.' f:'.sprintf("%.4f", $_->frequency).' c:'.$_->count}
    sort {$a->population->name cmp $b->population->name || $a->allele cmp $b->allele}
    @$al5
  ],
  [
    'p:ESP6500:AA a:- f:0.9740 c:4153',
    'p:ESP6500:AA a:T f:0.0260 c:111',
    'p:ESP6500:EA a:- f:0.9687 c:7982',
    'p:ESP6500:EA a:T f:0.0313 c:258'
  ],
  'get ESP allele frequency from VCF'
);

my $cdba = $multi->get_DBAdaptor('core');
my $sa = $cdba->get_SliceAdaptor();
my $slice = $sa->fetch_by_region('chromosome', 17); 
my $vf = $variation->get_all_VariationFeatures->[0];
my ($vcf_vf) = grep {$_->start == 3418000} @{$coll->get_all_VariationFeatures_by_Slice($slice)};
is(ref($vcf_vf->vcf_record), 'Bio::EnsEMBL::IO::Parser::VCF4Tabix', 'is Bio::EnsEMBL::IO::Parser::VCF4Tabix object');
my $al6 = $coll->get_all_Alleles_by_VariationFeature($vf, undef, $vcf_vf->vcf_record);

is_deeply(
  [
    map {'p:'.$_->population->name.' a:'.$_->allele.' f:'.sprintf("%.4f", $_->frequency).' c:'.$_->count}
    sort {$a->population->name cmp $b->population->name || $a->allele cmp $b->allele}
    @$al6
  ],
  [
    'p:ESP6500:AA a:- f:0.9740 c:4153',
    'p:ESP6500:AA a:T f:0.0260 c:111',
    'p:ESP6500:EA a:- f:0.9687 c:7982',
    'p:ESP6500:EA a:T f:0.0313 c:258'
  ],
  'get ESP allele frequency from VCF use Bio::EnsEMBL::Variation::VCFVariationFeature as argument'
);


done_testing();

