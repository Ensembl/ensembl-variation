#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2026] EMBL-European Bioinformatics Institute
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

use FileHandle;
use Bio::EnsEMBL::IO::Parser::VCF4Tabix;
use Bio::EnsEMBL::IO::Parser::VCF4;
use Bio::EnsEMBL::Registry;
use Data::Dumper;

my $chrom = $ENV{'LSB_JOBINDEX'};

if ($chrom == 23) {
  $chrom = 'X';
}
if ($chrom == 24) {
  $chrom = 'Y';
}

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
  -host => 'ensembldb.ensembl.org',
  -user => 'anonymous',
  -DB_VERSION => 95
);

my $vdba = $registry->get_DBAdaptor('human', 'variation');
$vdba->dbc->reconnect_when_lost(1);
my $cdba = $registry->get_DBAdaptor('human', 'core');
$cdba->dbc->reconnect_when_lost(1);

my $vfa = $vdba->get_VariationFeatureAdaptor;
my $slice_adaptor = $cdba->get_SliceAdaptor;

my $dir = '/variation/data/gnomAD/v2.1/grch37/genomes/';

my $fh = FileHandle->new("/remap_gnomad/genomes/unique_mappings/chrom$chrom.txt", 'r');

my $fh_results =  FileHandle->new("/remap_gnomad/genomes/unique_mappings/$chrom.vcf", 'w');

while (<$fh>) {
  chomp;
  my ($chrom38, $seq_region_id, $start38, $end38, $allele_string38, $variation_name, $variation_id, $chrom37, $start37, $end37, $allele_string37) = split/\t/;
  my ($ref, $alt) = split('/', $allele_string37);
  my ($ids, $qual, $filter, $info) = @{get_37_vcf_entry($chrom37, $start37, $end37, $ref, $alt)};
  ($ref, $alt) = split('/', $allele_string38);
  my ($vcf_pos, $vcf_ref, $vcf_alt) = @{get_vcf_coords_alleles($chrom38, $start38, $end38, $ref, $alt)};
  #CHROM POS ID REF ALT QUAL FILTER INFO
  if (!$vcf_pos) {
    print STDERR $_, "\n";
  }
  print $fh_results join("\t", $chrom38, $vcf_pos, $ids, $vcf_ref, $vcf_alt, $qual, $filter, $info), "\n";
}

$fh->close();
$fh_results->close();

sub get_37_vcf_entry {
  my ($chrom, $start, $end, $ref, $alt) = @_;
  ($start, $end) = ($start - 1, $end + 1) if ($start > $end);
  my $vcf_file = "$dir/gnomad.genomes.r2.1.sites.chr$chrom\_noVEP.vcf.gz";
  my $parser = Bio::EnsEMBL::IO::Parser::VCF4Tabix->open($vcf_file);
  $parser->seek($chrom, $start, $end);
  while ($parser->next) {
    my ($ref, $alts, $info) = (
      $parser->get_reference,
      $parser->get_alternatives,
      $parser->get_raw_info,
    );

    my ($chr, $start, $end) = (
      $parser->get_seqname,
      $parser->get_raw_start,
      $parser->get_raw_start,
    );
    $end += length($ref) - 1;
    my $is_indel = 0;
    foreach my $alt_allele(@$alts) {
      $is_indel = 1 if length($alt_allele) != length($ref);
    }
    if(scalar @$alts > 1) {
      if($is_indel) {

        # find out if all the alts start with the same base
        # ignore "*"-types
        my %first_bases = map {substr($_, 0, 1) => 1} grep {!/\*/} ($ref, @$alts);

        if(scalar keys %first_bases == 1) {
          $ref = substr($ref, 1) || '-';
          $start++;

          my @new_alts;

          foreach my $alt_allele(@$alts) {
            $alt_allele = substr($alt_allele, 1) unless $alt_allele =~ /\*/;
            $alt_allele = '-' if $alt_allele eq '';
            push @new_alts, $alt_allele;
          }

          $alts = \@new_alts;
        }
      }
    }
    elsif($is_indel) {
      # insertion or deletion (VCF 4+)
      if(substr($ref, 0, 1) eq substr($alts->[0], 0, 1)) {

        # chop off first base
        $ref = substr($ref, 1) || '-';
        $alts->[0] = substr($alts->[0], 1) || '-';

        $start++;
      }
    }

    my $vcf_allele_string = "$ref/".join('/', @$alts);
    if ($vcf_allele_string eq "$ref/$alt") {
      # get all VCF fields #CHROM POS ID REF ALT QUAL FILTER INFO
      my $ids = $parser->get_raw_IDs;
      my $qual = $parser->get_raw_score;
      my $filter = $parser->get_raw_filter_results;
      $parser->close;
      return [$ids, $qual, $filter, $info];
    }
  }
  $parser->close;
  return [];
}

sub get_vcf_coords_alleles {
  my ($chrom, $start, $end, $ref, $alt) = @_;
  my $slice = $slice_adaptor->fetch_by_region('chromosome', $chrom);
  my $vf = Bio::EnsEMBL::Variation::VariationFeature->new_fast({
    start          => $start,
    end            => $end,
    allele_string  => "$ref/$alt",
    strand         => 1,
    map_weight     => 1,
    adaptor        => $vfa,
    variation_name => 'variation_name',
    chr            => $chrom,
    slice          => $slice,
  });
  my $line = $vf->to_VCF_record();
  my ($vcf_pos, $vcf_ref, $vcf_alt) = ($line->[1], $line->[3], $line->[4]);
  return [$vcf_pos, $vcf_ref, $vcf_alt];
}


