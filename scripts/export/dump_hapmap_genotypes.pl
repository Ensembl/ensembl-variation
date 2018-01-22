#!/usr/bin/env perl
=head1 LICENSE
  Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
  Copyright [2016-2018] EMBL-European Bioinformatics Institute
  Licensed under the Apache License, Version 2.0 (the "License");
  you may not use this file except in compliance with the License.
  You may obtain a copy of the License at
       http://www.apache.org/licenses/LICENSE-2.0
  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.
=cut
=head1 CONTACT
 Please email comments or questions to the public Ensembl
 developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.
 Questions may also be sent to the Ensembl help desk at
 <http://www.ensembl.org/Help/Contact>.
=cut

use strict;
use warnings;

# fetch all hapmap populations and individuals
# get all genotypes
# dump variation, genotypes and (population frequency?)

use Bio::EnsEMBL::Registry;
use FileHandle;
use Getopt::Long;
use List::MoreUtils qw(first_index);

my $config = {};

GetOptions(
  $config,
  'registry=s',
  'species=s',
  'dir=s',
  'frequency_dir=s',
  'human_gvf_file=s',
  'population_gvf_file=s',
  'variation_ids_dir=s',
  'in_vcf_file=s',
  'out_vcf_file=s',
  'genotype_cache_dir=s',
  'assembly=s',
  'mode=s',
  'host=s',
  'user=s',
  'port=s',
) or die "Error: Failed to parse command line arguments\n";

=begin
Dump all genotypes for hapmap populations to a VCF file
Use populations:
CSHL-HAPMAP:HAPMAP-ASW 90
CSHL-HAPMAP:HAPMAP-CHB 90
CSHL-HAPMAP:HAPMAP-CHD 100
CSHL-HAPMAP:HAPMAP-GIH 100
CSHL-HAPMAP:HAPMAP-LWK 100 
CSHL-HAPMAP:HAPMAP-MEX 90 
CSHL-HAPMAP:HAPMAP-MKK 180
CSHL-HAPMAP:HAPMAP-TSI 100
CSHL-HAPMAP:HapMap-CEU 185
CSHL-HAPMAP:HapMap-HCB 48
CSHL-HAPMAP:HapMap-JPT 93
CSHL-HAPMAP:HapMap-YRI 185

submitter_handle: CSHL-HAPMAP

Steps:
dump_subsnp_ids
  - dump distinct variation_id and subsnp_id from allele where frequency_submitter_handle=2 (CSHL-HAPMAP)
  - 5037297 variation_ids/subsnp_ids
dump_genotypes
  - use all sample ids from the set of hapmap populations
  - dump all genotypes from compressed_genotype_var if subsnp_id is linked to CSHL-HAPMAP and sample_id in set
    of samples from hapmap populations
  - dump per chrom
  - dump looks like this: variation_id subsnp_id sample_id genotype_code 
correct_ychrom_genotypes
- fix problem with genotypes on Y chrom
dump_allele_frequencies
  - dump population allele frequencies per hapmap population (ASW, CHB, ..)
  - dump data from the allele table if frequency_submitter_handle=2  
  - this is dumped and added to the final vcf file info field
dump_gvf
  - use Homo_sapiens.gvf which contains variation_ids in attributes column
  - create new gvf file HAPMAP.gvf containing only hapmap variants retrieved from previous steps
  - add population allele frequencies to HAPMAP.gvf, only report frequency for alt alleles
19      344053  rs4897853       T       C       .       .       dbSNP_147;TSA=SNV;E_Freq;E_Hapmap;E_1000G;ASW-ALT=C;ASW-REF=T;ASW-FREQ=0.897959;CHB-ALT=C;CHB-REF=T;CHB-FREQ=0.634146;CHD-ALT=C;CHD-REF=T;CHD-FREQ=0.682353;GIH-ALT=C;GIH-REF=T;GIH-FREQ=0.840909;LWK-ALT=C;LWK-REF=T;LWK-FREQ=0.955556;MEX-ALT=C;MEX-REF=T;MEX-FREQ=0.83;MKK-ALT=C;MKK-REF=T;MKK-FREQ=0.891608;TSI-ALT=C;TSI-REF=T;TSI-FREQ=0.727273;CEU-ALT=C;CEU-REF=T;CEU-FREQ=0.761062;HCB-ALT=C;HCB-REF=T;HCB-FREQ=0.651163;JPT-ALT=C;JPT-REF=T;JPT-FREQ=0.773256;YRI-ALT=C;YRI-REF=T;YRI-FREQ=0.902655;Reference_seq=T;Variant_seq=C;variation_id=3279699;allele_string=T,C;AA=C 
gvf2vcf use script: ensembl-variation/scripts/misc/release/gvf2vcf.pl
  - parse HAPMAP.gvf to HAPMAP.vcf
add_genotypes_to_vcf
  - add genotypes to vcf file
=end
=cut

my $registry = 'Bio::EnsEMBL::Registry';

if ($config->{registry}) {
  $registry->load_all($config->{registry});
} else {
  $registry->load_registry_from_db(
    -host => $config->{host},
    -port => $config->{port},
    -user => $config->{user},
  );
}

my $species = $config->{species} || 'human';
my $dbc = $registry->get_DBAdaptor($species, 'variation')->dbc;

print_subsnp_ids('CSHL-HAPMAP') if ($config->{mode} eq 'print_subsnp_ids');
dump_genotypes('CSHL-HAPMAP') if ($config->{mode} eq 'dump_genotypes');
correct_ychrom_genotypes() if ($config->{mode} eq 'correct_ychrom_genotypes');
dump_allele_frequencies('CSHL-HAPMAP') if ($config->{mode} eq 'dump_allele_frequencies');
dump_gvf() if ($config->{mode} eq 'dump_gvf');
add_genotypes_to_vcf() if ($config->{mode} eq 'add_genotypes_to_vcf');
sub print_subsnp_ids {
  my $handle = shift;
  my $fh = FileHandle->new($config->{dir} . 'Subsnp_ids.txt', 'w');

  my $sth = $dbc->prepare(qq{
    SELECT distinct a.variation_id, a.subsnp_id
    FROM allele a, submitter_handle sh
    WHERE sh.handle_id = a.frequency_submitter_handle
    AND sh.handle = ?
  }, {mysql_use_result => 1});

  $sth->execute($handle);
  my ($variation_id, $subsnp_id);
  $sth->bind_columns(\$variation_id, \$subsnp_id);
  while ($sth->fetch) {
    print $fh "$variation_id\t$subsnp_id\n";
  }
  $sth->finish;
  $fh->close;
}

sub dump_genotypes {
  my $handle = shift;

  my $subsnp_ids = {};
  my $fh = FileHandle->new($config->{dir}. "/Subsnp_ids.txt", 'r'); 
  while (<$fh>) {
    chomp;
    my ($variation_id, $subsnp_id) = split/\t/;
    $subsnp_ids->{$subsnp_id} = 1;
  }
  $fh->close();
  my $slice_adaptor = $registry->get_adaptor($species, 'core', 'slice');
  my $seq_regions = $slice_adaptor->fetch_all('chromosome');
  my $population_adaptor = $registry->get_adaptor($species, 'variation', 'population');
  my $populations = {
    'CSHL-HAPMAP:HAPMAP-ASW' => 'ASW',
    'CSHL-HAPMAP:HAPMAP-CHB' => 'CHB',
    'CSHL-HAPMAP:HAPMAP-CHD' => 'CHD',
    'CSHL-HAPMAP:HAPMAP-GIH' => 'GIH',
    'CSHL-HAPMAP:HAPMAP-LWK' => 'LWK',
    'CSHL-HAPMAP:HAPMAP-MEX' => 'MEX',
    'CSHL-HAPMAP:HAPMAP-MKK' => 'MKK',
    'CSHL-HAPMAP:HAPMAP-TSI' => 'TSI',
    'CSHL-HAPMAP:HapMap-CEU' => 'CEU',
    'CSHL-HAPMAP:HapMap-HCB' => 'HCB',
    'CSHL-HAPMAP:HapMap-JPT' => 'JPT',
    'CSHL-HAPMAP:HapMap-YRI' => 'YRI',
  };

  my $sample_ids = {};
  foreach my $population_name (keys %$populations) {
    my $population = $population_adaptor->fetch_by_name($population_name);
    my @samples = @{$population->get_all_Samples};
    print STDERR $population_name, " ", scalar @samples, "\n";
    foreach my $sample (@samples) {
      $sample_ids->{$sample->dbID} = 1;
    }
  }

  my $sth = $dbc->prepare(qq{
    SELECT c.variation_id, c.subsnp_id, c.genotypes
    FROM seq_region sr, variation_feature vf, compressed_genotype_var c
    WHERE c.variation_id = vf.variation_id
    AND vf.seq_region_id = sr.seq_region_id
    AND sr.name = ?
  }, {mysql_use_result => 1});

  my $dir = $config->{dir};
  foreach my $seq_region (@$seq_regions) {
    my $seq_region_name = $seq_region->seq_region_name;
    my $fh = FileHandle->new($config->{dir} . "/genotypes/$seq_region_name.txt", 'w');
    $sth->execute($seq_region_name);

    my ($variation_id, $subsnp_id, $genotypes);
    $sth->bind_columns(\($variation_id, $subsnp_id, $genotypes));

    while ($sth->fetch) {
      next if (!$subsnp_ids->{$subsnp_id});
      my @genotypes = unpack("(ww)*", $genotypes);
      my @gt_codes = ();
      while (@genotypes) {
        my $sample_id = shift @genotypes;
        my $gt_code = shift @genotypes;
        if ($sample_ids->{$sample_id}) {
          print $fh join("\t", $variation_id, $subsnp_id, $sample_id, $gt_code), "\n";
        }
      }
    }
    $fh->close();
  }
}

sub correct_ychrom_genotypes {
  my $ychrom = $config->{dir} . "/genotypes/Y.txt";
  my $prevYchrom = $config->{dir} . "/genotypes/prevY.txt";
  `mv $ychrom $prevYchrom`;

  my $sample_ids = {};
  my $gtc_adaptor = $registry->get_adaptor($species, 'variation', 'genotypecode');
  my $gtcodes = {};

  my $fh = FileHandle->new($prevYchrom, 'r');
  while (<$fh>) {
    chomp;
    my ($variation_id, $subsnp_id, $sample_id, $gt_code_id) = split/\t/;
    $gtcodes->{$gt_code_id} = 1;
  }
  $fh->close;

  my $map_gtcodes = {};
  my $old2new = {};

  foreach my $gtcode_id (keys %$gtcodes) {
    my $gtc = $gtc_adaptor->fetch_by_dbID($gtcode_id);
    my @gt_alleles = @{$gtc->genotype};
    if ($gt_alleles[0] eq $gt_alleles[1]) {
      $map_gtcodes->{$gt_alleles[0]} = $gtcode_id;
    }
  }

  $gtcodes = $gtc_adaptor->fetch_all;

  foreach my $gtcode (@$gtcodes) {
    my @gt_alleles = @{$gtcode->genotype};
    if (scalar @gt_alleles == 1) {
      if ($map_gtcodes->{$gt_alleles[0]}) {
        $old2new->{$map_gtcodes->{$gt_alleles[0]}} = $gtcode->dbID;
      }
    }
  }

  my $fh_new = FileHandle->new($ychrom, 'w');

  $fh = FileHandle->new($prevYchrom, 'r');

  while (<$fh>) {
    chomp;
    my ($variation_id, $subsnp_id, $sample_id, $gt_code_id) = split/\t/;
    if ($old2new->{$gt_code_id}) {
      my $new_code_id = $old2new->{$gt_code_id};
      print $fh_new "$variation_id\t$subsnp_id\t$sample_id\t$new_code_id\n";
    }
  }
  $fh->close;
  $fh_new->close;

}

sub dump_allele_frequencies {
  my $handle = shift;
  my $pa = $registry->get_adaptor($species, 'variation', 'population');
  my $populations = {
    'CSHL-HAPMAP:HAPMAP-ASW' => 'ASW',
    'CSHL-HAPMAP:HAPMAP-CHB' => 'CHB',
    'CSHL-HAPMAP:HAPMAP-CHD' => 'CHD',
    'CSHL-HAPMAP:HAPMAP-GIH' => 'GIH',
    'CSHL-HAPMAP:HAPMAP-LWK' => 'LWK',
    'CSHL-HAPMAP:HAPMAP-MEX' => 'MEX',
    'CSHL-HAPMAP:HAPMAP-MKK' => 'MKK',
    'CSHL-HAPMAP:HAPMAP-TSI' => 'TSI',
    'CSHL-HAPMAP:HapMap-CEU' => 'CEU',
    'CSHL-HAPMAP:HapMap-HCB' => 'HCB',
    'CSHL-HAPMAP:HapMap-JPT' => 'JPT',
    'CSHL-HAPMAP:HapMap-YRI' => 'YRI',
  };

  my $fhs = {};

  foreach my $name (keys %$populations) {
    my $short_name = $populations->{$name};
    my $population = $pa->fetch_by_name($name);
    my $dbid = $population->dbID();
    my $fh = FileHandle->new( $config->{dir} . "/allele_frequencies/$short_name.txt", 'w');
    $fhs->{$dbid} = $fh;
  }

  my $sth = $dbc->prepare(qq{
    SELECT a.variation_id, ac.allele, a.frequency, a.population_id
    FROM allele a, allele_code ac, submitter_handle sh
    WHERE a.allele_code_id = ac.allele_code_id
    AND a.frequency_submitter_handle = sh.handle_id
    AND sh.handle = ?
  }, {mysql_use_result => 1});
  $sth->execute($handle);
  my ($variation_id, $allele, $frequency, $population_id);
  $sth->bind_columns(\($variation_id, $allele, $frequency, $population_id));
  while ($sth->fetch) {
    if ($population_id) {
      if ($fhs->{$population_id}) {
        my $fh = $fhs->{$population_id};
        print $fh join("\t", ($variation_id, $allele, $frequency)), "\n";
      }
    }
  }

  foreach my $dbid (keys %$fhs) {
    my $fh = $fhs->{$dbid};
    $fh->close;
  }

}

sub add_genotypes_to_vcf {
  my $cache = {};
  my $in_vcf_file = $config->{in_vcf_file};
  my $out_vcf_file = $config->{out_vcf_file};
  my $genotype_cache_dir = $config->{genotype_cache_dir};

  my $population_adaptor = $registry->get_adaptor($species, 'variation', 'population');
  my $populations = {
    'CSHL-HAPMAP:HAPMAP-ASW' => 'ASW',
    'CSHL-HAPMAP:HAPMAP-CHB' => 'CHB',
    'CSHL-HAPMAP:HAPMAP-CHD' => 'CHD',
    'CSHL-HAPMAP:HAPMAP-GIH' => 'GIH',
    'CSHL-HAPMAP:HAPMAP-LWK' => 'LWK',
    'CSHL-HAPMAP:HAPMAP-MEX' => 'MEX',
    'CSHL-HAPMAP:HAPMAP-MKK' => 'MKK',
    'CSHL-HAPMAP:HAPMAP-TSI' => 'TSI',
    'CSHL-HAPMAP:HapMap-CEU' => 'CEU',
    'CSHL-HAPMAP:HapMap-HCB' => 'HCB',
    'CSHL-HAPMAP:HapMap-JPT' => 'JPT',
    'CSHL-HAPMAP:HapMap-YRI' => 'YRI',
  };

  my $sample_ids = {};
  my $sample_id_2_name = {};
  my $sample_name_2_id = {};

  foreach my $population_name (keys %$populations) {
    my $population = $population_adaptor->fetch_by_name($population_name);
    my @samples = @{$population->get_all_Samples};
    foreach my $sample (@samples) {
      $sample_ids->{$sample->dbID} = 1;
      $sample_id_2_name->{$sample->dbID} = $sample->name;
      $sample_name_2_id->{$sample->name} = $sample->dbID;
    }
  }
  my @sample_names = sort keys %$sample_name_2_id;
  my $gtc_adaptor = $registry->get_adaptor($species, 'variation', 'genotypecode');
  my $gtcs = {};

  my $prev_chrom = -1;
  my $fh_in = FileHandle->new($in_vcf_file, 'r');

  my $fh_out = FileHandle->new($out_vcf_file, 'w');

  print $fh_out join("\n",
    '##fileformat=VCFv4.1',
    '##fileDate=20170104',
    '##reference=GRCh38',
    '##source=Genotypes for samples in HAPMAP populations: CSHL-HAPMAP:HAPMAP-ASW,CSHL-HAPMAP:HAPMAP-CHB,CSHL-HAPMAP:HAPMAP-CHD,CSHL-HAPMAP:HAPMAP-GIH,CSHL-HAPMAP:HAPMAP-LWK,CSHL-HAPMAP:HAPMAP-MEX,CSHL-HAPMAP:HAPMAP-MKK,CSHL-HAPMAP:HAPMAP-TSI,CSHL-HAPMAP:HapMap-CEU,CSHL-HAPMAP:HapMap-HCB,CSHL-HAPMAP:HapMap-JPT,CSHL-HAPMAP:HapMap-YRI',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
  ), "\n";

  print $fh_out join("\t", (qw/#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT/, @sample_names)), "\n";

  while (<$fh_in>) {
    chomp;
    next if (/^#/);
    my @values = split("\t", $_);
    my $chrom = $values[0];
    my $rs = $values[2];
    my $ref = $values[3];
    my @alts = split(',', $values[4]);
    my $info = get_info($values[7]);

    my $variation_id = $info->{variation_id};
    my $allele_string = $info->{allele_string};
    my @alleles = split(',', $allele_string);
    if ($chrom ne $prev_chrom) {
      print STDERR $chrom, "\n";
      $cache = {};
      $cache = update_genotype_cache("$genotype_cache_dir/$chrom.txt");
      print STDERR scalar keys %$cache, "\n";
    }
    $prev_chrom = $chrom; 

    my $sample_gts = $cache->{$variation_id};
    my @vcf_gts = ();
    foreach my $sample (@sample_names) {
      my $sample_id = $sample_name_2_id->{$sample};
      my $gt_code = $sample_gts->{$sample_id};

      if ($gt_code) {
        my $gtc = $gtcs->{$gt_code};
        if (!$gtc) {
          $gtc = $gtc_adaptor->fetch_by_dbID($gt_code);
          $gtcs->{$gt_code} = $gtc;
        }
        my $sep = ($gtc->phased) ? '|' : '/';
        my @gt_alleles = @{$gtc->genotype};
        my @indices = ();
        foreach my $allele (@gt_alleles) {
          my $index = first_index { $_ eq $allele } @alleles;
          push @indices, $index;
        }
        push @vcf_gts, join($sep, @indices);
      } else {
        if ($chrom eq 'Y') {
         push @vcf_gts, '.'; # different for X and Y
        } else {
          push @vcf_gts, './.'; # different for X and Y
        }
      }
    }

    print $fh_out join("\t", $values[0], $values[1], $values[2], $values[3], $values[4], $values[5], $values[6], '', 'GT', @vcf_gts), "\n";

  }

  $fh_in->close;
  $fh_out->close;
}

sub get_info {
  my $info = shift;
  my $hash = {};
  foreach (split/;/, $info) {
    if (/=/) {
      my ($key, $value) = split/=/;  
      $hash->{$key} = $value;
    }
  }
  return $hash;
}

sub update_genotype_cache {
  my $cache_file = shift;
  my $fh = FileHandle->new($cache_file, 'r');
  my $cache = {};
  while (<$fh>) {
    chomp;
    my ($variation_id, $subsnp_id, $sample_id, $gt_code) = split/\t/;
    $cache->{$variation_id}->{$sample_id} = $gt_code;
  }
  $fh->close();
  return $cache;
}

sub dump_gvf {
  my $cache = {};

  my $human_gvf_file = $config->{human_gvf_file};
  my $population_gvf_file = $config->{population_gvf_file};
  my $variation_ids_dir = $config->{variation_ids_dir};

  my $prev_seq_id = -1;
  my $fh_gvf = FileHandle->new($human_gvf_file, 'r');
  my $fh_out = FileHandle->new($population_gvf_file, 'w');
  while (<$fh_gvf>) {
    chomp;
    if (/^#/) {
      print $fh_out $_, "\n";
    } else {
      my $line = $_;
      my $gvf_line = get_gvf_line($line);
      my $seq_id = $gvf_line->{seq_id};
      unless ("$seq_id" eq "$prev_seq_id") {
        $cache = population_frequencies("$variation_ids_dir/$seq_id.txt");
        $prev_seq_id = $seq_id;
      }
      my $allele_string = $gvf_line->{attributes}->{allele_string};
      my $variation_id  = $gvf_line->{attributes}->{variation_id};
      next unless ($cache->{$variation_id});

      my @vf_alleles = split(',', $allele_string);

      if (my $population_allele_hash = $cache->{$variation_id}) {
        for my $short_name (keys %$population_allele_hash) {
          my $allele_hash = $population_allele_hash->{$short_name};
          my @alleles = ();
          my @freqs = ();

          for my $allele (keys %$allele_hash) {
            my $freq = $allele_hash->{$allele};
            if ($allele eq $vf_alleles[0]) {
              $gvf_line->{attributes}->{"$short_name-REF"} = $allele;
              if ($freq == 1) {
                if ((scalar keys %$allele_hash) == 1) {
                  push @alleles, $vf_alleles[1];
                  push @freqs, '0.0';
                }
              }
            } else {
              push @alleles, $allele;
              push @freqs, $freq;
            }
          }
          $gvf_line->{attributes}->{"$short_name-ALT"} = join(',', @alleles);
          if (grep {$_ =~ /\d/} @freqs) {
            $gvf_line->{attributes}->{"$short_name-FREQ"} = join(',', @freqs);
          }
        }
    }
    delete $gvf_line->{attributes}->{global_minor_allele_frequency};
    $line = join("\t", map {$gvf_line->{$_}} (
      'seq_id',
      'source',
      'type',
      'start',
      'end',
      'score',
      'strand',
      'phase'));
      my $attributes = join(";", map{"$_=$gvf_line->{attributes}->{$_}"} keys %{$gvf_line->{attributes}});
      print $fh_out $line, "\t", $attributes, "\n";
    }
  } # while
  $fh_gvf->close();
  $fh_out->close();
} 

sub population_frequencies {
  my $seq_id_file = shift;

  my $fh_var_ids = FileHandle->new($seq_id_file, 'r');
  my $variation_ids = {};
  while (<$fh_var_ids>) {
    chomp;
    my ($variation_id, $rest) = split("\t", $_, 2);
    $variation_ids->{$variation_id} = 1;
  }
  $fh_var_ids->close();

  my $frequency_dir = $config->{frequency_dir};
  my $cache = {};

  foreach my $short_name (qw/ASW CHB CHD GIH LWK MEX MKK TSI CEU HCB JPT YRI/) {
    my $fh = FileHandle->new("$frequency_dir/$short_name.txt", 'r');
    while (<$fh>) {
      chomp;
      my ($variation_id, $allele, $frequency) = split(/\t/);
      if ($variation_ids->{$variation_id}) {
        $frequency ||= 0;
        $cache->{$variation_id}->{$short_name}->{$allele} = $frequency;
      }
    }
    $fh->close();
  }
  return $cache;
}

sub get_gvf_line {
  my $line = shift;
  my $gvf_line = {};
  my @header_names = qw/seq_id source type start end score strand phase/;
  my @header_values = split(/\t/, $line);
  my $attrib = pop @header_values;

  for my $i (0 .. $#header_names) {
    $gvf_line->{$header_names[$i]} = $header_values[$i];
  }

  my @attributes = split(';', $attrib);
  foreach my $attribute (@attributes) {
    my ($key, $value) = split('=', $attribute);
    if ($value) {
      $gvf_line->{attributes}->{$key} = $value;
    }
  }
  return $gvf_line;
}
