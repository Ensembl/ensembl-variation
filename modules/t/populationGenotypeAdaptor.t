# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2023] EMBL-European Bioinformatics Institute
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
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Test::MultiTestDB;
use FileHandle;
my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');
my $vdba = $multi->get_DBAdaptor('variation');

#my $reg = 'Bio::EnsEMBL::Registry';
#$reg->load_all("$Bin/testnew.ensembl.registry");
#warn "using $Bin/testnew.ensembl.registry\n";
#my $vdba = $reg->get_DBAdaptor('human', 'variation');

my $pgta = $vdba->get_PopulationGenotypeAdaptor();

my $va = $vdba->get_VariationAdaptor();
my $pa = $vdba->get_PopulationAdaptor();


my $pgt = $pgta->fetch_by_dbID(12);

is($pgt->variation->name, 'rs1162', 'variation name');
is($pgt->allele(1), 'A', 'allele 1');
is($pgt->allele(2), 'A', 'allele 2');
is($pgt->frequency(), 0.4, 'frequency');
is($pgt->population->name, 'KWOK:HydatidiformMoles', 'population name');


my $variation_name = 'rs144235347';
#my $variation_name = 'rs7698608';
my $variation = $va->fetch_by_name($variation_name);

my $pgts = $pgta->fetch_all_by_Variation($variation);
is(scalar @$pgts, 22, 'number of pgts for variation');

my $hash;
foreach (@$pgts) {
    push @{$hash->{$_->population->name}}, $_->allele(1) . "|" . $_->allele(2);
}

## on the fly checks via compressed_genotype_var
is(join(',', sort @{$hash->{'1000GENOMES:phase_1_GBR'}}), 'C|C', 'pgts for GBR');
is(join(',', sort @{$hash->{'1000GENOMES:phase_1_LWK'}}), 'C|C,C|T', 'pgts for LWK');

my $population = $pa->fetch_by_name('PCJB:CUDAS');
my @pgts = sort {$a->dbID() <=> $b->dbID()} @{$pgta->fetch_all_by_Population($population)};

is(scalar @pgts, 3, 'number of pgts for population');
is($pgts[0]->variation()->name(), 'rs2255888', 'variation name');
is($pgts[0]->allele(1), 'T', 'allele1');
is($pgts[0]->allele(2), 'C', 'allele2');
is($pgts[0]->frequency(), 0.352252, 'frequency');
is($pgts[0]->population()->name(), 'PCJB:CUDAS', 'population name');

# store method
my $pg_hash = {
  genotype   => ['A','T'],
  frequency  => 0.12345,
  population => $pa->fetch_by_name('1000GENOMES:phase_1_GBR'),
  subsnp     => 12345,
  count      => 12345,
  variation  => $variation
};

my $pg1 = Bio::EnsEMBL::Variation::PopulationGenotype->new_fast($pg_hash);

ok($pgta->store($pg1), "store");

$pg_hash->{population} = $pa->fetch_by_name('1000GENOMES:phase_1_LWK');
my $pg2 = Bio::EnsEMBL::Variation::PopulationGenotype->new_fast($pg_hash);

my $hash2;
%$hash2 = %$pg_hash;
$hash2->{frequency} = 1 - $hash2->{frequency};
my $pg3 = Bio::EnsEMBL::Variation::PopulationGenotype->new_fast($hash2);

ok($pgta->store_multiple([$pg2, $pg3]), "store_multiple");

my $fetched = $variation->get_all_PopulationGenotypes();

ok((grep {$_->frequency && $_->frequency == 0.12345} @$fetched), "fetch stored pg 1");
ok((grep {$_->frequency && $_->frequency == 0.12345} grep {$_->population && $_->population->name eq '1000GENOMES:phase_1_LWK'} @$fetched), "fetch stored pg 2");
ok((grep {$_->frequency && $_->frequency == 0.87655} grep {$_->population && $_->population->name eq '1000GENOMES:phase_1_LWK'} @$fetched), "fetch stored pg 3");

my $fh = FileHandle->new;
my $tmpfile = "$Bin\/$$\_population_genotype.txt";
$fh->open(">$tmpfile") or die("ERROR: Could not write to $tmpfile\n");
ok($pgta->store_to_file_handle($pg1, $fh), "store_to_file_handle");
$fh->close();

open FH, $tmpfile;
my @lines = <FH>;
close FH;

unlink($tmpfile);

ok(scalar @lines == 1 && $lines[0] =~ /12345/, "check dumped data");

# read ESP population genotype frequencies from VCF
my $dir = $multi->curr_dir();
ok($vdba->vcf_config_file($dir.'/vcf_config.json') eq $dir.'/vcf_config.json', "DBAdaptor vcf_config_file");
my $vca = $vdba->get_VCFCollectionAdaptor();
my $coll = $vca->fetch_by_id('esp_GRCh37');
my $temp = $coll->filename_template();
$temp =~ s/###t\-root###/$dir/;
$coll->filename_template($temp);
$pgta->db->use_vcf(1);

$variation = $va->fetch_by_name('rs35099512');
my $pg4 = $pgta->fetch_all_by_Variation($variation);

is_deeply(
  [
    map {'p:'.$_->population->name.' gt:'.$_->genotype_string.' f:'.sprintf("%.4f", $_->frequency).' c:'.$_->count}
    sort {$a->population->name cmp $b->population->name || $a->genotype_string cmp $b->genotype_string}
    @$pg4
  ],
  [
    'p:ESP6500:AA gt:-|- f:0.9498 c:2025',
    'p:ESP6500:AA gt:-|T f:0.0483 c:103',
    'p:ESP6500:AA gt:T|T f:0.0019 c:4',
    'p:ESP6500:EA gt:-|- f:0.9408 c:3876',
    'p:ESP6500:EA gt:-|T f:0.0558 c:230',
    'p:ESP6500:EA gt:T|T f:0.0034 c:14'
  ],
  'get ESP population genotype frequency from VCF'
);

my ($esp_ea_population) = grep {$_->name eq 'ESP6500:EA'} @{$coll->get_all_Populations};
ok($esp_ea_population && $esp_ea_population->isa('Bio::EnsEMBL::Variation::Population'), "grep by population name from get_all_Populations");

my $pg5 = $pgta->fetch_all_by_Variation($variation, $esp_ea_population);
is_deeply(
  [
    map {'p:'.$_->population->name.' gt:'.$_->genotype_string.' f:'.sprintf("%.4f", $_->frequency).' c:'.$_->count}
    sort {$a->population->name cmp $b->population->name || $a->genotype_string cmp $b->genotype_string}
    @$pg5
  ],
  [
    'p:ESP6500:EA gt:-|- f:0.9408 c:3876',
    'p:ESP6500:EA gt:-|T f:0.0558 c:230',
    'p:ESP6500:EA gt:T|T f:0.0034 c:14'
  ],
  'get ESP population genotype frequency from VCF for ESP6500:EA'
);

my $cdba = $multi->get_DBAdaptor('core');
my $sa = $cdba->get_SliceAdaptor();
my $slice = $sa->fetch_by_region('chromosome', 17);

my $vf = Bio::EnsEMBL::Variation::VariationFeature->new(
  -seq_region_name => 17,
  -start => 3417997,
  -end   => 3417996,
  -allele_string => 'CACCA/ACCA',
  -slice => $slice,
  -strand => 1,
);
my ($vcf_vf) = grep {$_->start == 3417997} @{$coll->get_all_VariationFeatures_by_Slice($slice)};
is(ref($vcf_vf->vcf_record), 'Bio::EnsEMBL::IO::Parser::VCF4Tabix', 'is Bio::EnsEMBL::IO::Parser::VCF4Tabix object');
my $pg6 = $coll->get_all_PopulationGenotypes_by_VariationFeature($vf, undef, $vcf_vf->vcf_record);

is_deeply(
  [
    map {'p:'.$_->population->name.' gt:'.$_->genotype_string.' f:'.sprintf("%.4f", $_->frequency).' c:'.$_->count}
    sort {$a->population->name cmp $b->population->name || $a->genotype_string cmp $b->genotype_string}
    @$pg6
  ],
  [
    'p:ESP6500:AA gt:ACCA|ACCA f:0.0014 c:3',
    'p:ESP6500:AA gt:ACCA|CACCA f:0.0455 c:97',
    'p:ESP6500:AA gt:CACCA|CACCA f:0.9507 c:2026',
    'p:ESP6500:AA gt:CACCA|T f:0.0023 c:5',
    'p:ESP6500:EA gt:ACCA|ACCA f:0.0024 c:10',
    'p:ESP6500:EA gt:ACCA|CACCA f:0.0514 c:212',
    'p:ESP6500:EA gt:ACCA|T f:0.0002 c:1',
    'p:ESP6500:EA gt:CACCA|CACCA f:0.9423 c:3885',
    'p:ESP6500:EA gt:CACCA|T f:0.0024 c:10',
    'p:ESP6500:EA gt:T|T f:0.0012 c:5',
  ],
  'get ESP population genotype frequency from VCF use Bio::EnsEMBL::Variation::VCFVariationFeature as argument'
);

$coll = $vca->fetch_by_id('esp_GRCh37_with_error');
$temp = $coll->filename_template();
$temp =~ s/###t\-root###/$dir/;
$coll->filename_template($temp);
$vf = $variation->get_all_VariationFeatures->[0];
warns_like { $coll->get_all_PopulationGenotypes_by_VariationFeature($vf); } qr/Unexpected genotype variable. Expected values are R, A1, A2, An at/, 'Warn on unsupported genotype string in VCF file.';

done_testing();
