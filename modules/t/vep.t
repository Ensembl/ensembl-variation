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
use Test::Exception;
use Data::Dumper;
use FindBin qw($Bin);
use File::Path qw(remove_tree);;

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
  get_version_data
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

# get version data
my $vd = get_version_data($base_config);
is_deeply($vd, {
  'sift' => 'sift5.2.2',
  'polyphen' => '2.2.2',
  'COSMIC' => '71',
  'ESP' => '20140509',
  'gencode' => 'GENCODE',
  'HGMD-PUBLIC' => '20142',
  'genebuild' => '2014-07',
  'regbuild' => '13.0',
  'assembly' => 'GRCh38.p2',
  'dbSNP' => '138',
  'ClinVar' => '201410'
}, "cache version data");

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
    'STRAND' => -1,
    'IMPACT' => 'MODERATE',
    'FLAGS' => 'cds_end_NF',
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

is_deeply($exp, $cons->[2], "get_all_consequences 2");

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
  format     => 'ensembl',
  allele_number => 1,
  gencode_basic => 1,
  variant_class => 1,
});

$cons = get_all_consequences($config, [$vf]);

ok($cons && scalar @$cons == 3, "get_all_consequences - everything 1");

$exp = {
  'IMPACT' => 'MODERATE',
  'SYMBOL' => 'MRPL39',
  'DOMAINS' => 'hmmpanther:PTHR11451,Gene3D:3.10.20.30,Superfamily_domains:SSF81271',
  'SYMBOL_SOURCE' => 'HGNC',
  'ENSP' => 'ENSP00000404426',
  'ALLELE_NUM' => 1,
  'PolyPhen' => 'probably_damaging(0.975)',
  'BIOTYPE' => 'protein_coding',
  'UNIPARC' => 'UPI0000E5A387',
  'AA_MAF' => 'C:0',
  'SIFT' => 'deleterious(0)',
  'STRAND' => -1,
  'HGNC_ID' => 'HGNC:14027',
  'HGVSc' => 'ENST00000419219.1:c.275N>G',
  'HGVSp' => 'ENSP00000404426.1:p.Ala92Gly',
  'TREMBL' => 'C9JG87',
  'EA_MAF' => 'C:0.0001',
  'VARIANT_CLASS' => 'SNV',
  'EXON' => '2/8',
  'TSL' => '5',
  'FLAGS' => 'cds_end_NF',
};

is_deeply($cons->[2]->{Extra}, $exp, "get_all_consequences - everything 2");

# regulatory
$config = copy_config($base_config, {
  regulatory => 1,
  cell_type  => ['HUVEC'],
  biotype    => 1,
});
($vf) = @{parse_line($config, '21 25562380 25562380 G/T +')};
$cons = get_all_consequences($config, [$vf]);

ok((grep {$_->{Extra}->{BIOTYPE} && $_->{Extra}->{BIOTYPE} eq 'promoter'} @$cons), "regulatory - type");
ok((grep {$_->{Extra}->{MOTIF_SCORE_CHANGE} && $_->{Extra}->{MOTIF_SCORE_CHANGE} == -0.017} @$cons), "regulatory - motif score");

# regulatory SV
$config = copy_config($base_config, { regulatory => 1 });
($vf) = @{parse_line($config, '21 25560805 25568206 DEL +')};
$cons = get_all_consequences($config, [$vf]);
ok((grep {$_->{Feature} eq 'ENSR00000612089'} @$cons), "regulatory - SV overlap");

## input formats

# ensembl SV
$config = copy_config($base_config);
($vf) = @{parse_line($config, '21 25587759 25587769 DEL + del')};
$cons = get_all_consequences($config, [$vf]);
ok((grep {$_->{Consequence} =~ /feature_truncation/} @$cons), "ensembl format - SV del");

$config = copy_config($base_config);
($vf) = @{parse_line($config, '21 25587759 25587769 DUP + del')};
$cons = get_all_consequences($config, [$vf]);
ok((grep {$_->{Consequence} =~ /feature_elongation/} @$cons), "ensembl format - SV dup");

# intergenic SV
$config = copy_config($base_config);
($vf) = @{parse_line($config, '20 25587759 25587769 DEL + del')};
$cons = get_all_consequences($config, [$vf]);
ok((grep {$_->{Consequence} =~ /intergenic/} @$cons), "ensembl format - intergenic SV");

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
21      25587759        sv_del          .       <DEL>   .       .       SVTYPE=DEL;END=25587769;CIPOS=5,5;CIEND=5,5
21      25587759        sv_dup          .       <DUP>   .       .       SVTYPE=DUP;END=25587769
1 1 nochr G A . . .};

my @lines = split("\n", $input);

$config = copy_config($base_config);
($vf) = @{parse_line($config, $lines[0])};
ok($vf && $vf->allele_string eq 'C/T', "vcf format - parse_line");

# SV input
my @vfs = grep {validate_vf($config, $_)} map {@{parse_line($config, $_)}} @lines;
ok((grep {$_->isa('Bio::EnsEMBL::Variation::StructuralVariationFeature')} @vfs) == 2, "vcf format - SVs");

$config->{pick_allele} = 1;
$cons = get_all_consequences($config, \@vfs);

my @sv_cons = grep {$_->{Uploaded_variation} eq 'sv_del'} @$cons;
ok((grep {$_->{Consequence} =~ /feature_truncation/} @sv_cons), "vcf format - SV del cons");

@sv_cons = grep {$_->{Uploaded_variation} eq 'sv_dup'} @$cons;
ok((grep {$_->{Consequence} =~ /feature_elongation/} @sv_cons), "vcf format - SV dup cons");

# vcf deletion
$config = copy_config($base_config, {allow_non_variant => 1, vcf => 1});
($vf) = @{parse_line($config, '21 25606454 test GC G')};
ok($vf && $vf->allele_string eq 'C/-' && $vf->start == 25606455 && $vf->end == $vf->start, "vcf format - deletion");

# vcf insertion
$config = copy_config($base_config, {allow_non_variant => 1, vcf => 1});
($vf) = @{parse_line($config, '21 25606454 test G GC')};
ok($vf && $vf->allele_string eq '-/C' && $vf->start == 25606455 && $vf->end == 25606454, "vcf format - insertion");

# vcf multiple alleles
$config = copy_config($base_config);
($vf) = @{parse_line($config, '21 25606454 test G C,T')};
ok($vf && $vf->allele_string eq 'G/C/T', "vcf format - multiple alleles");
$cons = get_all_consequences($config, [$vf]);
ok($cons && scalar @$cons == 6, "vcf format - multiple alleles cons");

# vcf mixed allele types
$config = copy_config($base_config);
($vf) = @{parse_line($config, '21 25606454 test G C,TT')};
ok($vf && $vf->allele_string eq 'G/C/TT', "vcf format - mixed allele types 1");

$config = copy_config($base_config);
($vf) = @{parse_line($config, '21 25606454 test G GC,GT')};
ok($vf && $vf->allele_string eq '-/C/T', "vcf format - mixed allele types 2");

# vcf non variant
$config = copy_config($base_config, {allow_non_variant => 1, vcf => 1});
($vf) = @{parse_line($config, '21 25606454 test G . . . CSQ=A')};
$cons = get_all_consequences($config, [$vf]);
ok($cons && ${$cons->[0]} =~ /21\s+25606454\s+test\s+G\s+./, "vcf format - non variant");

# check csq removed
ok(${$cons->[0]} !~ /CSQ\=A/, "vcf format - existing CSQ removed");

# vcf weird "*" type
$config = copy_config($base_config);
($vf) = @{parse_line($config, '21 25606454 test G C,* . . .')};
is($vf->allele_string, 'G/C/*', 'vcf format - *-type 1');

$config = copy_config($base_config);
($vf) = @{parse_line($config, '21 25606454 test GC G,* . . .')};
is($vf->allele_string, 'C/-/*', 'vcf format - *-type 1 indel');

# this type is produced by GATK, we have to force VEP to parse it as VCF as it's not in the standard
$config = copy_config($base_config, {format => 'vcf'});
($vf) = @{parse_line($config, '21 25606454 test G C,<DEL:*> . . .')};
is($vf->allele_string, 'G/C/<DEL:*>', 'vcf format - *-type 2');


# use minimal to reduce allele strings
($vf) = @{parse_line($config, '21 25606454 test GA GT')};
is($vf->allele_string, 'GA/GT', "minimal - dont use");

$config = copy_config($base_config, {minimal => 1});
($vf) = @{parse_line($config, '21 25606454 test GA GT')};
is($vf->allele_string, 'A/T', "minimal - use");
is($vf->start, 25606455, "minimal - start coord adjusted");
is($vf->end, 25606455, "minimal - end coord adjusted");

($vf) = @{parse_line($config, '21 25606454 test GAC GTC')};
is($vf->allele_string, 'A/T', "minimal - trim both ends");

($vf) = @{parse_line($config, '21 25741665 test CAGAAGAAAG TAGAAGAAAG,C')};
my %by_allele = map {$_->{Allele} => $_} grep {$_->{Feature} eq 'ENST00000400075'} @{get_all_consequences($config, [$vf])};
is($by_allele{'T'}->{Consequence}, 'missense_variant', "minimal - complex (SNP)");
is($by_allele{'-'}->{Consequence}, 'inframe_deletion,splice_region_variant', "minimal - complex (del)");
is($by_allele{'-'}->{Extra}->{MINIMISED}, 1, "minimal - extra flag");

$config = copy_config($base_config, {minimal => 1, allele_number => 1});
($vf) = @{parse_line($config, '21 25606454 test G GC,C')};
$cons = get_all_consequences($config, [$vf]);
ok(defined($cons->[0]->{Extra}->{ALLELE_NUM}), "minimal - allele num defined");

%by_allele = map {$_->{Consequence} => $_} grep {$_->{Feature} eq 'ENST00000419219'} @$cons;
is(scalar keys %{{map {$_->{Allele} => 1} values %by_allele}}, 1, "minimal - allele num where two alts resolve to same allele (check allele)");
is($by_allele{frameshift_variant}->{Extra}->{ALLELE_NUM}, 1, "minimal - allele num where two alts resolve to same allele (frameshift)");
is($by_allele{missense_variant}->{Extra}->{ALLELE_NUM}, 2, "minimal - allele num where two alts resolve to same allele (missense)");

$config = copy_config($base_config, {allele_number => 1});
($vf) = @{parse_line($config, '21 25606454 test G *,C')};
$cons = get_all_consequences($config, [$vf]);
%by_allele = map {$_->{Allele} => $_} grep {$_->{Feature} eq 'ENST00000419219'} @$cons;
is($by_allele{C}->{Extra}->{ALLELE_NUM}, 2, "allele_num - * first allele");

# vcf individual data
$config = copy_config($base_config, {
  individual => ['all'],
  ind_cols => {
    'A' => 9,
    'B' => 10,
  }
});
($vf) = @{parse_line($config, qq{21 25587758 rs116645811 G A . . . GT 1|1 0|0})};
$cons = get_all_consequences($config, [$vf]);
ok($cons && (grep {$_->{Extra}->{IND} eq 'A'} @$cons) && !(grep {$_->{Extra}->{IND} eq 'B'} @$cons), "vcf format - individual data");

# vcf process_ref_homs
$config = copy_config($base_config, {
  individual => ['all'],
  ind_cols => {
    'A' => 9,
    'B' => 10,
  },
  process_ref_homs => 1,
});
@vfs = @{parse_line($config, qq{21 25587758 rs116645811 G A . . . GT 1|1 0|0})};
ok(@vfs && $vfs[1]->{individual} eq 'B', "vcf format - individual data process ref homs");

# vcf GP
$config = copy_config($base_config, { gp => 1 });
($vf) = @{parse_line($config, qq{1 1 test G C . . GP=21:25606454})};
ok($vf && $vf->start == 25606454 && $vf->{chr} eq '21', "vcf format - gp");


# pileup
$config = copy_config($base_config);
($vf) = grep {validate_vf($config, $_)} @{parse_line($config, 'chr21 25606454 G C')};
ok($vf && $vf->allele_string eq 'G/C', "pileup format - parse_line");
$cons = get_all_consequences($config, [$vf]);
ok($cons && scalar @$cons == 3, "pileup format - consequences");

# pileup indel

# invalid format
$config = copy_config($base_config);
dies_ok { parse_line($config, 'a b') } "invalid format";


## validate_vf
$config = copy_config($base_config);
# delete($config->{quiet});   # uncomment this to see warnings
$vf = Bio::EnsEMBL::Variation::VariationFeature->new_fast({
  chr => 1,
  start => 123,
  end => 123,
  allele_string => 'A/C',
});

# invalid coord
$vf->{start} = '1234a567';
ok(!validate_vf($config, $vf), "validate_vf - invalid coord");

# invalid allele string
$vf->{start} = 123;
$vf->{allele_string} = '9';
ok(!validate_vf($config, $vf), "validate_vf - invalid allele_string");

# allele looks like insertion
$vf->{allele_string} = '-/A';
ok(!validate_vf($config, $vf), "validate_vf - alleles look like insertion");

# start > end + 1
$vf->{allele_string} = 'A/C';
$vf->{start} = 125;
ok(!validate_vf($config, $vf), "validate_vf - start > end + 1");

# alleles not compatible with coords
# $vf->{allele_string} = 'AA/C';
# $vf->{start} = 123;
# ok(!validate_vf($config, $vf), "validate_vf - alleles not compatible with coords");


## other options

# frequency filtering
$input = qq{21 25000248 25000248 C/G + test1
21 25000264 25000264 A/G + test2};

$config = copy_config($base_config, {
  check_existing => 1,
  check_frequency => 1,
  freq_pop => '1kg_eas',
  freq_freq => 0.03,
  freq_gt_lt => 'lt',
  freq_filter => 'include',
});
$cons = get_all_consequences($config, [map {@{parse_line($config, $_)}} split("\n", $input)]);
my %ex = map {$_->{Uploaded_variation} => 1} @$cons;
ok(!$ex{test1} && $ex{test2}, "check frequency 1");
is($cons->[0]->{Extra}->{FREQS}, '1kg_eas:0.0010', "check frequency 2");

$config = copy_config($base_config, {
  check_existing => 1,
  check_frequency => 1,
  freq_pop => '1kg_eas',
  freq_freq => 0.03,
  freq_gt_lt => 'lt',
  freq_filter => 'exclude',
});
$cons = get_all_consequences($config, [map {@{parse_line($config, $_)}} split("\n", $input)]);
%ex = map {$_->{Uploaded_variation} => 1} @$cons;
ok($ex{test1} && !$ex{test2}, "check frequency 3");

$config = copy_config($base_config, {
  check_existing => 1,
  check_frequency => 1,
  freq_pop => '1kg_eas',
  freq_freq => 0.03,
  freq_gt_lt => 'gt',
  freq_filter => 'include',
});
$cons = get_all_consequences($config, [map {@{parse_line($config, $_)}} split("\n", $input)]);
%ex = map {$_->{Uploaded_variation} => 1} @$cons;
ok($ex{test1} && !$ex{test2}, "check frequency 4");

# summary
$config = copy_config($base_config, { summary => 1 });
($vf) = @{parse_line($config, '21 25587758 rs116645811 G A . . .')};
$cons = get_all_consequences($config, [$vf]);
my %got = map {$_ => 1} split(',', $cons->[0]->{Consequence});
%ex = (
  missense_variant => 1,
  intron_variant => 1,
  upstream_gene_variant => 1,
);
is_deeply(\%got, \%ex, "summary");

# most severe
$config = copy_config($base_config, { most_severe => 1 });
($vf) = @{parse_line($config, '21 25587758 rs116645811 G A . . .')};
$cons = get_all_consequences($config, [$vf]);
%got = map {$_ => 1} split(',', $cons->[0]->{Consequence});
%ex = (
  missense_variant => 1,
);
is_deeply(\%got, \%ex, "most_severe");

# most severe regulatory
$config = copy_config($base_config, {
  regulatory => 1,
});
($vf) = @{parse_line($config, '21 25487468 25487468 A/T +')};
$cons = get_all_consequences($config, [$vf]);
%got = map {$_ => 1} split(',', $cons->[0]->{Consequence});
%ex = (
  regulatory_region_variant => 1,
);
is_deeply(\%got, \%ex, "most_severe regulatory");

# flag pick
$config = copy_config($base_config, { flag_pick => 1 });
($vf) = @{parse_line($config, '21 25587758 rs116645811 G A . . .')};
$cons = get_all_consequences($config, [$vf]);
ok($cons && (grep {$_->{Extra}->{PICK}} grep {$_->{Consequence} eq 'missense_variant'} @$cons), "flag pick");

# flag pick allele
$config = copy_config($base_config, { flag_pick_allele => 1 });
($vf) = @{parse_line($config, '21 25587758 rs116645811 G A,C . . .')};
$cons = get_all_consequences($config, [$vf]);
ok($cons && (grep {$_->{Extra}->{PICK}} @$cons) == 2, "flag pick allele");

# fork
$input = qq{21      25607440        rs61735760      C       T       .       .       .
21      25606638        rs3989369       A       G       .       .       .
21      25606478        rs75377686      T       C       .       .       .
21      25603925        rs7278284       C       T       .       .       .};
@lines = split("\n", $input);
$config = copy_config($base_config, { fork => 2, pick => 1, vcf => 1 });
@vfs = grep {validate_vf($config, $_)} map {@{parse_line($config, $_)}} @lines;
$exp = [qw(25607440 25606638 25606478 25603925)];

$cons = get_all_consequences($config, \@vfs);
ok($cons && scalar @$cons == 4, "fork");

my $got = [map {(split(/\s+/, $$_))[1]} @$cons];
is_deeply($got, $exp, "fork - order preserved");

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
$config = copy_config($base_config, {
  json => 1,
  rest => 1,
});

($vf) = @{parse_line($config, '21 25606454 25606454 G/C +')};
$cons = get_all_consequences($config, [$vf]);
ok($cons && $cons->[0]->{most_severe_consequence} eq 'missense_variant', "json output");


$config = copy_config($base_config, {
  json => 1,
  rest => 1,
  hgvs => 1,
});

($vf) = @{parse_line($config, '21 25606453 25606453 G/C +')};
$cons = get_all_consequences($config, [$vf]);

ok($cons && $cons->[0]->{transcript_consequences}->[2]->{hgvsp} eq 'ENSP00000404426.1:p.Ala92=', "json HGVSp no escaping");


# check
$config = copy_config($base_config, {
  json => 1,
  rest => 1,
  check_existing => 1,
  pubmed => 1,
});
($vf) = @{parse_line($config, '21 25227941 25227941 T/C +')};
$cons = get_all_consequences($config, [$vf, $vf]);
is_deeply(
  $cons->[0]->{colocated_variants},
  [{
    'phenotype_or_disease' => 1,
    'eas_allele' => 'C',
    'amr_maf' => '0.0432',
    'strand' => 1,
    'id' => 'rs41504145',
    'sas_allele' => 'C',
    'sas_maf' => '0.0133',
    'allele_string' => 'T/C',
    'amr_allele' => 'C',
    'minor_allele_freq' => '0.127',
    'afr_allele' => 'C',
    'eas_maf' => '0.003',
    'afr_maf' => '0.4289',
    'eur_maf' => '0.0229',
    'end' => 25227941,
    'eur_allele' => 'C',
    'minor_allele' => 'C',
    'start' => 25227941,
    'pubmed' => [17903302],
  }],
  'json output - duplicate position'
);

# solr xml
$config = copy_config($base_config, {
  solr => 1,
  fields => ['Consequence', 'SIFT', 'GMAF'],
  pick_allele => 1,
  sift => 'b',
  check_existing => 1,
  gmaf => 1,
});

($vf) = @{parse_line($config, '21 25606454 25606454 G/C +')};
$cons = get_all_consequences($config, [$vf]);
ok($cons && ${$cons->[0]} =~ /\<field name\=\"Consequence\">missense_variant\<\/field\>/, "solr xml output");


## CUSTOM FILES
if(`which tabix` =~ /tabix/) {
  $input = qq{21      25592860        rs10576 T       C       .       .       .
21      25592836        rs1135638       G       A       .       .       .
21      25587758        rs116645811     G       A       .       .       .};
  $config = copy_config($base_config, {
    custom => [
      {
        file   => "$Bin\/testdata/test.bed.gz",
        name   => 'testbed',
        type   => 'overlap',
        format => 'bed',
        coords => 0
      },
      {
        file   => "$Bin\/testdata/test.gff.gz",
        name   => 'testgff',
        type   => 'overlap',
        format => 'gff',
        coords => 0
      },
      {
        file   => "$Bin\/testdata/test.vcf.gz",
        name   => 'testvcf',
        type   => 'exact',
        format => 'vcf',
        coords => 0,
        fields => ['ATTR'],
      }
    ],
    pick => 1,
  });
  @vfs = map {@{parse_line($config, $_)}} split("\n", $input);
  $cons = get_all_consequences($config, \@vfs);

  my %by_var = map {$_->{Uploaded_variation} => $_->{Extra}} @$cons;
  ok($by_var{rs116645811}->{testbed} eq 'bed1' && $by_var{rs1135638}->{testbed} eq 'bed2', "custom - bed");

  ok($by_var{rs116645811}->{testgff} eq 'gtf1', "custom - gff");

  ok(
    $by_var{rs116645811}->{testvcf} eq 'vcf1' &&
    $by_var{rs1135638}->{testvcf} eq 'vcf2' &&
    $by_var{rs10576}->{testvcf} eq 'vcf3',
    "custom - vcf"
  );
}
else {
  print STDERR "# tabix not found, skipping custom file tests\n";
}


## DATABASE
###########

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new('homo_sapiens');
my $vdb = $multi->get_DBAdaptor('variation');
my $cdb = $multi->get_DBAdaptor('core');
my $rdb = $multi->get_DBAdaptor('funcgen');

# make DB config
$config = copy_config($base_config, {
  
  database => 1,
  hgvs => 1,
  
  # core adaptors
  sa  => $cdb->get_SliceAdaptor,
  ta  => $cdb->get_TranscriptAdaptor,
  ga  => $cdb->get_GeneAdaptor,
  csa => $cdb->get_CoordSystemAdaptor,
  
  # var adaptors
  va    => $vdb->get_VariationAdaptor,
  vfa   => $vdb->get_VariationFeatureAdaptor,
  tva   => $vdb->get_TranscriptVariationAdaptor,
  svfa  => $vdb->get_StructuralVariationFeatureAdaptor,
  pfpma => $vdb->get_ProteinFunctionPredictionMatrixAdaptor,
  
  # reg adaptors
  RegulatoryFeature_adaptor => $rdb->get_RegulatoryFeatureAdaptor,
  MotifFeature_adaptor      => $rdb->get_MotifFeatureAdaptor,
});
delete $config->{cache};

# parse HGVS
($vf) = @{parse_line($config, "2:g.46739212C>G")};
ok($vf && $vf->allele_string eq 'C/G', "parse_line HGVS");

$cons = get_all_consequences($config, [$vf]);
ok((grep {$_->{Extra}->{HGVSc} eq 'ENST00000522587.1:c.639G>C'} @$cons), "DB output 1");

# parse ID
delete $config->{format};
($vf) = @{parse_line($config, "rs2299222")};
ok($vf && $vf->variation_name eq 'rs2299222', "parse_line ID");

# check ref
delete $config->{format};
$config->{check_ref} = 1;
($vf) = @{parse_line($config, "2 46739212 46739212 C/G +")};
$vf->{slice} = get_slice($config, $vf->{chr}, undef, 1);

ok(validate_vf($config, $vf), "validate_vf - check ref pass");

$vf->{allele_string} = 'T/G';
ok(!validate_vf($config, $vf), "validate_vf - check ref fail");

delete $config->{check_ref};

# check svs
$config->{check_svs} = 1;
($vf) = @{parse_line($config, "8 7803895 7803895 C/T +")};
$vf->{slice} = get_slice($config, $vf->{chr}, undef, 1);
$cons = get_all_consequences($config, [$vf]);
ok($cons && $cons->[0]->{Extra}->{SV} eq 'esv89107', "check svs");

# map from non-toplevel coord
$config->{cache} = 1;
($vf) = grep {validate_vf($config, $_)} @{parse_line($config, "AC018682.4 1 1 C/G +")};
ok($vf && $vf->{chr} eq '2' && $vf->{original_chr} eq 'AC018682.4', "non-toplevel transform");

# failed map
($vf) = grep {validate_vf($config, $_)} @{parse_line($config, "mangledAC018682.4 1 1 C/G +")};
ok(!$vf, "non-toplevel transform fail");
delete $config->{cache};


# map to LRG

# regulatory
delete($config->{format});
$config->{regulatory} = 1;
($vf) = @{parse_line($config, "7 151409212 151409212 C/T +")};
$vf->{slice} = get_slice($config, $vf->{chr}, undef, 1);
$cons = get_all_consequences($config, [$vf]);

my ($rf_con) = grep {$_->{Feature_type} && $_->{Feature_type} eq 'RegulatoryFeature'} @$cons;
my ($mf_con) = grep {$_->{Feature_type} && $_->{Feature_type} eq 'MotifFeature'} @$cons;

is($rf_con->{Extra}->{BIOTYPE}, 'regulatory_region', "db - regulatory biotype");
is($rf_con->{Feature}, 'ENSR00000636355', "db - regulatory ID");
is($mf_con->{Feature}, 'PB0043.1', "db - motif ID");
is_deeply(
  $mf_con->{Extra},
  {
    'STRAND' => -1,
    'MOTIF_POS' => 16,
    'MOTIF_NAME' => 'Jaspar_Matrix_Max:PB0043.1',
    'HIGH_INF_POS' => 'N',
    'MOTIF_SCORE_CHANGE' => '0.010',
    'IMPACT' => 'MODIFIER',
  },
  "db - motif extra"
);

# build
$config = copy_config($config, {
  reg         => 'Bio::EnsEMBL::Registry',
  build       => 22,
  build_parts => 'tvr',
  strip       => 1,
  write_cache => 1,
  symbol      => 1,
  cell_type   => [1],
  freq_vcf    => [
    {
      pops => ['AFR','AMR'],
      prefixed_pops => ['AFR','AMR'],
      file => "$Bin\/testdata/freqs.vcf.gz",
    },
  ],
});

# this needs resetting as it might be updated in a newer version
delete $config->{cache_variation_cols};

$config->{dir} = "$Bin\/testdata/$$\_vep_cache/".$config->{species}."/".$config->{cache_version}."_".$config->{assembly};

build_full_cache($config);

$config = copy_config($base_config, {
  check_existing => 1,
  maf_1kg        => 1,
});
$config->{dir} =~ s/$Bin\/testdata\/vep-cache/$Bin\/testdata\/$$\_vep_cache/;

ok(read_cache_info($config), "read built cache info");

($vf) = grep {validate_vf($config, $_)} @{parse_line($config, "22 20876358 20876358 C/T +")};
$cons = get_all_consequences($config, [$vf]);
ok($cons && scalar @$cons == 1 && $cons->[0]->{Feature} eq 'ENST00000420225', "build - basic test");

ok($cons->[0]->{Extra}->{AFR_MAF} eq 'T:0.03' && $cons->[0]->{Extra}->{AMR_MAF} eq 'T:0.05', "build - freqs from vcf 1");
ok($cons->[0]->{Extra}->{CLIN_SIG} eq 'pathogenic', "build - clin_sig");

($vf) = grep {validate_vf($config, $_)} @{parse_line($config, "22 20876359 20876359 C/G +")};
$cons = get_all_consequences($config, [$vf]);
ok($cons->[0]->{Extra}->{AFR_MAF} eq 'T:0.03' && $cons->[0]->{Extra}->{AMR_MAF} eq 'T:0.05', "build - freqs from vcf 2");

($vf) = grep {validate_vf($config, $_)} @{parse_line($config, "22 20876360 20876360 T/G +")};
$cons = get_all_consequences($config, [$vf]);
ok($cons->[0]->{Extra}->{AFR_MAF} eq 'G:0.03' && $cons->[0]->{Extra}->{AMR_MAF} eq 'G:0.05', "build - freqs from vcf 3");

($vf) = @{parse_line($config, "22 50655788 50655788 G/C +")};
$cons = get_all_consequences(copy_config($config, {pick => 1}), [$vf]);
ok($cons->[0]->{Amino_acids} eq 'U/S', "build - selenocysteine seqEdit applied");

$config->{regulatory} = 1;
($vf) = grep {validate_vf($config, $_)} @{parse_line($config, "22 20001112 20001112 T/G +")};
$cons = get_all_consequences($config, [$vf]);
($rf_con) = grep {$_->{Feature_type} && $_->{Feature_type} eq 'RegulatoryFeature'} @$cons;
($mf_con) = grep {$_->{Feature_type} && $_->{Feature_type} eq 'MotifFeature'} @$cons;

# is($rf_con->{Extra}->{BIOTYPE}, 'regulatory_region', "build - regulatory biotype");
is($rf_con->{Feature}, 'ENSR00000672895', "build - regulatory ID");
is($mf_con->{Feature}, 'MA0139.1', "build - motif ID");
is_deeply(
  $mf_con->{Extra},
  {
    'STRAND' => -1,
    'MOTIF_POS' => 18,
    'MOTIF_NAME' => 'Jaspar_Matrix_CTCF:MA0139.1',
    'HIGH_INF_POS' => 'N',
    'MOTIF_SCORE_CHANGE' => '0.000',
    'IMPACT' => 'MODIFIER',
  },
  "build - motif extra"
);

# remove built cache
remove_tree("$Bin\/testdata/$$\_vep_cache");

done_testing();

sub copy_config {
  my $config = shift;
  my $copy   = shift;
  $copy ||= {};
  $copy->{$_} = $config->{$_} for keys %$config;
  return $copy;
}
