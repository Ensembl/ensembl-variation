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
use Data::Dumper;
use FindBin qw($Bin);

use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Variation::Utils::VEP qw(
  parse_line
  vf_to_consequences
  validate_vf
  convert_to_vcf
  get_all_consequences
  get_slice
  build_full_cache
  read_cache_info
  get_time
  debug
);

# configure
my $base_config = {};

open CONF, "$Bin\/vep.conf" or die "ERROR: Could not read from conf file $Bin\/test.conf\n";
while(<CONF>) {
  chomp;
  my ($k, $v) = split("\t", $_);
  $v =~ s/###t\-root###/$Bin/g;
  $base_config->{$k} = $v;
}
close CONF;

# read_cache_info
ok(read_cache_info($base_config), "read_cache_info");

# parse line
my $config = copy_config($base_config);
my ($vf) = @{parse_line($config, '21 25606454 25606454 G/C +')};
ok($vf && $vf->isa('Bio::EnsEMBL::Variation::VariationFeature'), "parse_line 1");
ok($vf->allele_string eq 'G/C', "parse_line 2");
ok($vf->class_SO_term eq 'SNV', "parse_line 3");

# validate_vf
ok(validate_vf($config, $vf), "validate_vf");

# get_all_consequences
my $cons = get_all_consequences($config, [$vf]);

ok($cons && scalar @$cons == 3, "get_all_consequences 1");

my $exp = {
  'Consequence' => 'missense_variant',
  'Extra' => {
    'STRAND' => -1
  },
  'Feature_type' => 'Transcript',
  'Uploaded_variation' => undef,
  'Existing_variation' => '-',
  'Allele' => 'C',
  'Gene' => 'ENSG00000154719',
  'CDS_position' => '275',
  'cDNA_position' => '284',
  'Protein_position' => '92',
  'Amino_acids' => 'A/G',
  'Feature' => 'ENST00000419219',
  'Codons' => 'gCc/gGc',
  'Location' => '21:25606454'
};

is_deeply($exp, $cons->[0], "get_all_consequences 2");

# make a copy of $config with loads switched on
$config = copy_config($base_config, {
  sift       => 'b',
  polyphen   => 'b',
  ccds       => 1,
  hgvs       => 1,
  symbol     => 1,
  numbers    => 1,
  domains    => 1,
  canonical  => 1,
  protein    => 1,
  biotype    => 1,
  gmaf       => 1,
  check_existing => 1,
  maf_1kg    => 1,
  maf_esp    => 1,
  pubmed     => 1,
  uniprot    => 1,
  tsl        => 1,
  format     => 'ensembl'
});

$cons = get_all_consequences($config, [$vf]);

ok($cons && scalar @$cons == 3, "get_all_consequences - everything 1");

$exp = {
  'SYMBOL' => 'MRPL39',
  'DOMAINS' => 'Superfamily_domains:SSF81271',
  'SYMBOL_SOURCE' => 'HGNC',
  'ENSP' => 'ENSP00000404426',
  'PolyPhen' => 'probably_damaging(0.975)',
  'BIOTYPE' => 'protein_coding',
  'UNIPARC' => 'UPI0000E5A387',
  'AA_MAF' => 'C:0',
  'SIFT' => 'deleterious(0.01)',
  'STRAND' => -1,
  'HGNC_ID' => 'HGNC:14027',
  'HGVSc' => 'ENST00000419219.1:c.275N>G',
  'HGVSp' => 'ENSP00000404426.1:p.Ala92Gly',
  'TREMBL' => 'C9JG87',
  'EA_MAF' => 'C:0.000116',
  'EXON' => '2/8',
  'TSL' => '5'
};

is_deeply($exp, $cons->[0]->{Extra}, "get_all_consequences - everything 2");

# regulatory
$config = copy_config($base_config, {
  regulatory => 1,
  cell_type  => ['HUVEC'],
  biotype    => 1,
});
($vf) = @{parse_line($config, '21 25487468 25487468 A/T +')};
$cons = get_all_consequences($config, [$vf]);

ok((grep {$_->{Extra}->{BIOTYPE} && $_->{Extra}->{BIOTYPE} eq 'promoter_flanking_region'} @$cons), "regulatory - type");
ok((grep {$_->{Extra}->{MOTIF_SCORE_CHANGE} && $_->{Extra}->{MOTIF_SCORE_CHANGE} == -0.022} @$cons), "regulatory - motif score");

# vcf
my $input = qq{21      25607440        rs61735760      C       T       .       .       .
21      25606638        rs3989369       A       G       .       .       .
21      25606478        rs75377686      T       C       .       .       .
21      25603925        rs7278284       C       T       .       .       .
21      25603910        rs7278168       C       T       .       .       .
21      25603832        rs116331755     A       G       .       .       .
21      25592893        rs1057885       T       C       .       .       .
21      25592860        rs10576         T       C       .       .       .
21      25592836        rs1135638       G       A       .       .       .
21      25587758        rs116645811     G       A       .       .       .
21      25587759        sv_del          .       <DEL>   .       .       SVTYPE=DEL;END=25587769
21      25587759        sv_dup          .       <DUP>   .       .       SVTYPE=DUP;END=25587769};

my @lines = split("\n", $input);

$config = copy_config($base_config);
($vf) = @{parse_line($config, $lines[0])};
ok($vf && $vf->allele_string eq 'C/T', "parse_line vcf");

my @vfs = grep {validate_vf($config, $_)} map {@{parse_line($config, $_)}} @lines;
ok(scalar @vfs == 12, "SVs pass validate_vf");

$config->{pick_allele} = 1;
$cons = get_all_consequences($config, \@vfs);

my @sv_cons = grep {$_->{Uploaded_variation} eq 'sv_del'} @$cons;
ok((grep {$_->{Consequence} =~ /feature_truncation/} @sv_cons), "SV del cons");

@sv_cons = grep {$_->{Uploaded_variation} eq 'sv_dup'} @$cons;
ok((grep {$_->{Consequence} =~ /feature_elongation/} @sv_cons), "SV dup cons");


## output formats

# vcf
$config = copy_config($base_config, {
  vcf    => 1,
  fields => ['Feature','Consequence'],
  per_gene => 1,
});

($vf) = @{parse_line($config, '21 25606454 25606454 G/C +')};
$cons = get_all_consequences($config, [$vf]);
ok($cons && ${$cons->[0]} =~ /CSQ\=ENST00000307301\|missense_variant/, "vcf output");

# gvf
$config = copy_config($base_config, {gvf => 1});

($vf) = @{parse_line($config, '21 25606454 25606454 G/C +')};
$vf->variation_name('test');
$cons = get_all_consequences($config, [$vf]);
ok($cons && ${$cons->[0]} =~ /missense_variant 0 mRNA ENST00000419219/, "gvf output");

# json
$config = copy_config($base_config, {json => 1, rest => 1});

($vf) = @{parse_line($config, '21 25606454 25606454 G/C +')};
$cons = get_all_consequences($config, [$vf]);
ok($cons && $cons->[0]->{most_severe_consequence} eq 'missense_variant', "json output");

# solr xml
$config = copy_config($base_config, {
  solr => 1,
  fields => ['Consequence'],
  summary => 1,
});

($vf) = @{parse_line($config, '21 25606454 25606454 G/C +')};
$cons = get_all_consequences($config, [$vf]);
ok($cons && ${$cons->[0]} =~ /\<field name\=\"Consequence\">missense_variant\<\/field\>/, "solr xml output");



## DATABASE
###########

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');
my $vdb = $multi->get_DBAdaptor('variation');
my $cdb = $multi->get_DBAdaptor('core');

# make DB config
$config = copy_config($base_config, {
  
  database => 1,
  hgvs => 1,
  
  # core adaptors
  sa => $cdb->get_SliceAdaptor,
  ta => $cdb->get_TranscriptAdaptor,
  ga => $cdb->get_GeneAdaptor,
  
  # var adaptors
  va  => $vdb->get_VariationAdaptor,
  vfa => $vdb->get_VariationFeatureAdaptor,
  tva => $vdb->get_TranscriptVariationAdaptor,
});
delete $config->{$_} for qw(cache dir);

# parse HGVS
($vf) = @{parse_line($config, "2:g.46739212C>G")};
ok($vf && $vf->allele_string eq 'C/G', "parse_line HGVS");

$cons = get_all_consequences($config, [$vf]);
ok((grep {$_->{Extra}->{HGVSc} eq 'ENST00000522587.1:c.639G>C'} @$cons), "DB output 1");

# parse ID
delete $config->{format};
($vf) = @{parse_line($config, "rs2299222")};
ok($vf && $vf->variation_name eq 'rs2299222', "parse_line ID");

done_testing();

sub copy_config {
  my $config = shift;
  my $copy   = shift;
  $copy ||= {};
  $copy->{$_} = $config->{$_} for keys %$config;
  return $copy;
}
