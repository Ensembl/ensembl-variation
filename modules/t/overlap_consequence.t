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
use File::Spec::Functions;
use List::MoreUtils qw(first_index uniq);

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Variation::Utils::Config;
use Bio::EnsEMBL::Variation::Utils::Constants;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');
my $cdba = $multi->get_DBAdaptor('core');
my $vdba = $multi->get_DBAdaptor('variation');

# check if all consequences have unique ranks
my @csqs = @{Bio::EnsEMBL::Variation::Utils::Config::OVERLAP_CONSEQUENCES};
my %ranks;
for (@csqs) {
  $ranks{$_->{rank}} = [] unless $ranks{$_->{rank}};
  push($ranks{$_->{rank}}, $_->{SO_term});
}

for my $rank ( %ranks ) {
  next unless $ranks{$rank};
  ok(scalar @{$ranks{$rank}} == 1,
     "Ranks must be unique: rank $rank assigned to consequence(s) " .
       join(", ", @{$ranks{$rank}}));
}

# check if consequences are identical between Constants.pm and Config.pm
# Constants.pm must be generated with ensembl-variation/scripts/misc/create_config_consts.pl
my %const_csqs = %{Bio::EnsEMBL::Variation::Utils::Constants::OVERLAP_CONSEQUENCES};
for my $csq (keys %const_csqs) {
  my @matches = grep { $_->{SO_term} eq $csq } @csqs;
  is_deeply($const_csqs{$csq}, $matches[0],
            "Identical consequences in Config.pm and Constants.pm");
}

# check if consequence is properly ranked based on ensembl-webcode's ConsequenceTable.pm
for my $match (grep { /ensembl-webcode/ } uniq @INC) {
  my $table = catdir($match, "EnsEMBL", "Web", "Document", "HTML", "ConsequenceTable.pm");
  if (-e $table) {
    open(my $fh, '<', $table) or next;

    # read all consequence terms in order from ConsequenceTable.pm
    my @table_terms;
    while(my $line = <$fh>){
      next unless $line =~ /'term' +=>/;
      my ($term) = $line =~ /=> '(.*)'/;
      push @table_terms, $term;
    }
    close($fh);

    # check if order is maintained
    for my $csq (@csqs) {
      my $index = first_index { $_ eq $csq->{SO_term} } @table_terms;
      ok($index > -1, "$csq->{SO_term} not listed in ConsequenceTable.pm");
      next unless $index > -1;
      cmp_ok($index + 1, "==", $csq->{rank},
             "$csq->{SO_term} rank must match order in ConsequenceTable.pm");
    }
  }
}

my $vfa = $vdba->get_VariationFeatureAdaptor;
my $slice_adaptor = $cdba->get_SliceAdaptor;
my $chrom = 13;
my $slice = $slice_adaptor->fetch_by_region('chromosome', $chrom);
my $vf = Bio::EnsEMBL::Variation::VariationFeature->new_fast({
    start          => 51483899,
    end            => 51484001,
    allele_string  =>  "GCCCTGCGGCGCGCGATGAGCACCTACCACTAGCGGAGCCGCGAGGGAGAGGCCGCGGCCCCTTCCCGTTGCCTGCGGCCACCGGCCGGCATTCAGAGCCCCT/-",
    strand         => 1,
    map_weight     => 1,
    slice          => $slice,
    adaptor        => $vfa,
});
my $consequence_hash = {};
my $vfos = $vf->get_all_VariationFeatureOverlaps;
foreach my $vfo (@{$vfos}) {
  my $vfoas = $vfo->get_all_alternate_VariationFeatureOverlapAlleles;
  foreach my $vfoa (@$vfoas) {
    my $consequences = $vfoa->get_all_OverlapConsequences;
    foreach my $consequence (@$consequences) {
      $consequence_hash->{$consequence->SO_term} = 1;
    }
  }
}
my $expected = '5_prime_UTR_variant,non_coding_transcript_exon_variant,regulatory_region_ablation,regulatory_region_variant';
my $got = join(',', sort keys %$consequence_hash);
ok($got eq $expected, "All overlap consequences for variation feature overlaps");

$vf = Bio::EnsEMBL::Variation::VariationFeature->new_fast({
    start          => 51483899,
    end            => 51484001,
    allele_string  =>  "GCCCTGCGGCGCGCGATGAGCACCTACCACTAGCGGAGCCGCGAGGGAGAGGCCGCGGCCCCTTCCCGTTGCCTGCGGCCACCGGCCGGCATTCAGAGCCCCT/-",
    strand         => 1,
    map_weight     => 1,
    slice          => $slice,
    adaptor        => $vfa,
});
$consequence_hash = {};
my $rfvs = $vf->get_all_RegulatoryFeatureVariations;
foreach my $rfv (@$rfvs) {
  my $rfvas = $rfv->get_all_alternate_RegulatoryFeatureVariationAlleles;
  foreach my $rfva (@$rfvas) {
    my $consequences = $rfva->get_all_OverlapConsequences;
    foreach my $consequence (@$consequences) {
      $consequence_hash->{$consequence->SO_term} = 1;
    }
  }
}

$expected = 'regulatory_region_ablation,regulatory_region_variant';
$got = join(',', sort keys %$consequence_hash);
ok($got eq $expected, "All overlap consequences for regulatory variation features");

done_testing();
