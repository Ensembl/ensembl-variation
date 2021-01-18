=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2021] EMBL-European Bioinformatics Institute

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

=head1 DESCRIPTION

This module contains functions used in the variant quality control process. 

=cut


package Bio::EnsEMBL::Variation::Utils::RemappingUtils;

use strict;
use warnings;

use base qw(Exporter);
use Bio::DB::Fasta;
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp expand);
use Bio::EnsEMBL::Variation::Utils::Sequence qw(revcomp_tandem);

our @EXPORT_OK = qw(map_variant filter_vf_mapping_first_round filter_vf_mapping_second_round qc_mapped_vf filter_svf_mapping filter_read_mapping);

my $NO_MAPPING = 'no mapping';
my $TOO_MANY_MAPPINGS = 'map weight > 5';
my $POOR_SCORE = 'does not pass alignment score';

my $failed_descriptions = {
  1 =>  'Variant maps to more than 3 different locations',
  2 =>  'None of the variant alleles match the reference allele',
  3 =>  'Variant has more than 3 different alleles',
  4 =>  'Loci with no observed variant alleles in dbSNP',
  5 =>  'Variant does not map to the genome',
  6 =>  'Variant has no genotypes',
  7 =>  'Genotype frequencies do not add up to 1',
  8 =>  'Variant has no associated sequence',
  9 =>  'Variant submission has been withdrawn by the 1000 genomes project due to high false positive rate',
  11 => 'Additional submitted allele data from dbSNP does not agree with the dbSNP refSNP alleles',
  12 => 'Variant has more than 3 different submitted alleles',
  13 => 'Alleles contain non-nucleotide characters',
  14 => 'Alleles contain ambiguity codes',
  15 => 'Mapped position is not compatible with reported alleles',
  16 => 'Flagged as suspect by dbSNP',
  17 => 'Variant can not be re-mapped to the current assembly',
  18 => 'Supporting evidence can not be re-mapped to the current assembly',
  19 => 'Variant maps to more than one genomic location',
  20 => 'Variant at first base in sequence',
};

my @copy_over_failure_reasons = (3, 4, 6, 7, 8, 9, 11, 12, 13, 14, 16, 18, 20);

sub map_variant {
  my $alignment         = shift;
  my $length_before_var = shift;
  my $length_var        = shift;

  my $indel_colocation = 0;
  my $flag_suspicious  = 0;
  my @cigar_string = @{$alignment->cigar_array};

  my $q_start = $alignment->query->start;
  my $q_end   = $alignment->query->end;
  my $t_start = $alignment->start;
  my $t_end   = $alignment->end;
  my ($snp_t_start, $snp_t_end);
  if ($q_start && $q_end && $t_start && $t_end) {
    $t_start = $t_start - $q_start + 1;
    $q_start = 1;  
  } else {
    $flag_suspicious = 1;    
    return [undef, undef, 0, $flag_suspicious];
  } 
  my $snp_q_pos = $length_before_var;     

  while (my $sub_string = shift @cigar_string) {
    my $operation = $sub_string->[0];
    my $count     = $sub_string->[1];
    my ($query_match_length, $target_match_length, $new_q_start, $new_t_start);
    if ($operation eq 'M' || $operation eq 'H' || $operation eq 'S') {
      $query_match_length  = $count;
      $target_match_length = $count;  
      $operation           = 'M';
    } elsif ($operation eq 'I') {
      $query_match_length  = $count;
      $target_match_length = 0;
    } elsif ($operation eq 'D') {   
      $query_match_length  = 0;
      $target_match_length = $count;
    }

    $new_q_start = $q_start + $query_match_length - 1 if ($operation eq 'M');
    $new_q_start = $q_start + $query_match_length if ($operation eq 'D' || $operation eq 'I');  

    $new_t_start = $t_start + $target_match_length - 1 if ($operation eq 'M');
    $new_t_start = $t_start + $target_match_length if ($operation eq 'D' || $operation eq 'I');


    if (($q_start <= $snp_q_pos and $snp_q_pos < $new_q_start) or 
        ($new_q_start < $snp_q_pos and $snp_q_pos <= $q_start)) {
      if ($operation eq 'M') {
        if ($length_var == 0) {
          $snp_t_start = $t_start + abs($snp_q_pos - $q_start);
          $snp_t_end = $snp_t_start + 1;
          ($snp_t_start, $snp_t_end) = ($snp_t_end, $snp_t_start);
        } else {
          $snp_t_start = $t_start + abs($snp_q_pos - $q_start + 1);
          $snp_t_end = $snp_t_start + $length_var - 1;
        } 
      } else { 
        if ($operation eq 'I' && ($new_q_start > $snp_q_pos)) {
          $indel_colocation = 1;
        }
        if ($length_var == 0) {
          $snp_t_start = $new_t_start;
          $snp_t_end = $snp_t_start + 1;
          ($snp_t_start, $snp_t_end) = ($snp_t_end, $snp_t_start);
        } else {
          $snp_t_start = $new_t_start + 1;
          $snp_t_end = $snp_t_start + $length_var - 1;
        }
      } 
    }
    if ($snp_t_start && $snp_t_end) {
      last;   
    }
    $q_start = $new_q_start;
    $t_start = $new_t_start;  
  } # end while
  return [$snp_t_start, $snp_t_end, $indel_colocation, $flag_suspicious];
}


sub filter_read_mapping {
  my $config = shift;
  my $feature_type_id = $config->{feature_type_id}; # e.g. 'phenotype_feature_id';
  my $file_init_feature = $config->{file_init_feature};
  my $file_filtered_mappings = $config->{file_filtered_mappings};
  my $file_mappings = $config->{file_mappings};

  my $fh_mappings = FileHandle->new($file_mappings, 'r');
  my $mapped = {};
  while (<$fh_mappings>) {
    chomp;
    #query_name seq_region_name start end strand map_weight edit_dist score_count algn_score
    #10137:patched    4       127441656       127441755       1       1       1.00000 241M
    #417:complete     4       127285626       127285726       1       293     0.99010 101M
    #10157:upstream   4       128525045       128525544       -1      1       0.99800 277M1I223M
    #10157:downstream 4       129826530       129827031       -1      1       0.99401 213M1D288M
    my ($query_name, $new_seq_name, $new_start, $new_end, $new_strand, $map_weight, $algn_score, $cigar) = split("\t", $_);
    my ($id, $type) = split(':', $query_name);
    $mapped->{$id}->{$type}->{$new_seq_name}->{"$new_start:$new_end:$new_strand"} = $algn_score;
  }
  $fh_mappings->close();

  my $fh_init_feature = FileHandle->new($file_init_feature, 'r');
  # this comes from the initial dump from phenotype_feature
  my $feature_data = {};
  while (<$fh_init_feature>) {
    chomp;
    my $data = read_line($_);
    $feature_data->{$data->{$feature_type_id}} = $data;
  }
  $fh_init_feature->close();
  my ($stats_failed, $stats_unique_map, $stats_multi_map);
  my $fh_filtered_mappings = FileHandle->new($file_filtered_mappings, 'w');
  foreach my $id (keys %$mapped) {
    my $init_data = $feature_data->{$id};
    my $old_chrom = $init_data->{seq_region_name};
    my $old_start = $init_data->{seq_region_start};
    my $old_end = $init_data->{seq_region_end};
    my $diff = $old_end - $old_start + 1;
    my $types = scalar keys %{$mapped->{$id}}; # upstream, downstream, patched, complete
    my $count_out = 0;
    my @output_lines = ();
    if ($types == 2) { # upstream and downstream mappings
      my $filtered = {};
      foreach my $seq_name (keys %{$mapped->{$id}->{upstream}}) {
        # find pairs
        my $filtered_upstream   = $mapped->{$id}->{upstream}->{$seq_name};
        my $filtered_downstream = $mapped->{$id}->{downstream}->{$seq_name};
        next unless ($filtered_upstream && $filtered_downstream);
        # match up upstream and downstream reads select pairs whose difference matches best with the size (diff) of the originial QTL
        $filtered = match_upstream_downstream($filtered_upstream, $filtered_downstream, $filtered, $diff, $seq_name);
      }
      # use prior knowledge from previous mapping and prefer same chromosome on which QLT was located in old assembly
      if ($filtered->{$old_chrom}) {
        filter_by_prior($filtered->{$old_chrom}, \@output_lines, {count_out => \$count_out, old_chrom => $old_chrom, id => $id});
      } else {
        # compute best score over all seq_names
        # Get best score for each seq_region
        # Get all seq_regions with best score
        # Get all mappings from seq_regions with best score
        filter_by_best_overall_score($filtered, \@output_lines, {count_out => \$count_out, old_chrom => $old_chrom, id => $id});
      }
    } elsif ($types == 1) { # complete or patched
      foreach my $type (keys %{$mapped->{$id}}) {
        my $mappings = get_patched_or_complete_mappings($mapped->{$id}->{$type}, $old_chrom);
        filter_by_score($mappings, \@output_lines, {count_out => \$count_out, id => $id});
      }
    } else {
      print STDERR "Error for id: $id. More than 2 types\n";
    }
    if ($count_out == 1) {
      $stats_unique_map++;
    } elsif ($count_out > 1 && $count_out <= 5) {
      $stats_multi_map++;
    } elsif ($count_out > 5) {
      $stats_failed++;
    } else {
      $stats_failed++;
    }
    if ($count_out >= 1 && $count_out <= 5) {
      foreach my $line (@output_lines) {
        print $fh_filtered_mappings $line;
      }
    }
  }
  $fh_filtered_mappings->close();
  return {
    unique_map => $stats_unique_map,
    multi_map => $stats_multi_map,
    failed => $stats_failed,
  };
}

sub filter_by_score {
  my $mappings = shift;
  my $output_lines = shift;
  my $hash = shift;
  my $count_out = $hash->{count_out};
  my $id = $hash->{id};
  my @sorted_mappings_by_score = sort {$mappings->{$a} <=> $mappings->{$b}} keys(%$mappings);
  my $threshold = $mappings->{$sorted_mappings_by_score[0]};
  foreach my $mapping (@sorted_mappings_by_score) {
    my $score = $mappings->{$mapping};
    if ($score <= $threshold) {
      $$count_out++;
      my ($chrom, $new_start, $new_end, $new_strand) = split(':', $mapping);
      push @$output_lines, "$id\t$chrom\t$new_start\t$new_end\t$new_strand\t$score\n";
    }
  }
}

sub get_patched_or_complete_mappings {
  my $mapped = shift;
  my $old_chrom = shift;
  my $mappings = {};
  if ($mapped->{$old_chrom}) {
    foreach my $coords (keys %{$mapped->{$old_chrom}}) {
      $mappings->{"$old_chrom:$coords"} = $mapped->{$old_chrom}->{$coords};
    }
  } else {
    foreach my $seq_name (keys %$mapped) {
      foreach my $coords (keys %{$mapped->{$seq_name}}) {
        $mappings->{"$seq_name:$coords"} = $mapped->{$seq_name}->{$coords};
      }
    }
  }
  return $mappings;
}

sub match_upstream_downstream {
  my $filtered_upstream = shift;
  my $filtered_downstream = shift;
  my $filtered = shift;
  my $diff = shift;
  my $seq_name = shift;
  foreach my $upstream_mapping (keys %$filtered_upstream) {
    my $upstream_score = $filtered_upstream->{$upstream_mapping};
    my ($up_start, $up_end, $up_strand) = split(':', $upstream_mapping);
    foreach my $downstream_mapping (keys %$filtered_downstream) {
      my $downstream_score = $filtered_downstream->{$downstream_mapping};
      my ($down_start, $down_end, $down_strand) = split(':', $downstream_mapping);
      my $mapping_diff = $down_end - $up_start + 1;
      my $diff_score = abs($mapping_diff - $diff);
      if ($up_strand == $down_strand) {
        if ($up_start > $down_end) {
          ($up_start, $down_end) = ($down_end, $up_start);
        }
        $filtered->{$seq_name}->{"$up_start:$down_end:$down_strand"} = $diff_score;
      }
    }
  }
  return $filtered;
}

sub filter_by_prior {
  my $mappings = shift;
  my $output_lines = shift;
  my $hash = shift;
  my $count_out = $hash->{count_out};
  my $old_chrom = $hash->{old_chrom};
  my $id = $hash->{id};
  my @sorted_diff_score = sort {$mappings->{$a} <=> $mappings->{$b}} keys(%$mappings); # score is defined as difference between previous length of QTL and new length of QTL, the smaller the better
  my $threshold = $mappings->{$sorted_diff_score[0]};
  foreach my $mapping (@sorted_diff_score) {
    my $score = $mappings->{$mapping};
    if ($score <= $threshold) { # consider all mappings with the same score
      $$count_out++;
      my ($new_start, $new_end, $new_strand) = split(':', $mapping);
      push @$output_lines, "$id\t$old_chrom\t$new_start\t$new_end\t$new_strand\t$score\n";
    }
  }
}

sub filter_by_best_overall_score {
  my $filtered = shift;
  my $output_lines = shift;
  my $hash = shift;
  my $count_out = $hash->{count_out};
  my $old_chrom = $hash->{old_chrom};
  my $id = $hash->{id};
  my $seq_region_2_best_score = {};
  # find best mapping score for each seq_region
  foreach my $seq_name (keys %$filtered) {
    my $mappings = $filtered->{$seq_name};
    my @sorted_mappings_by_score = sort {$mappings->{$a} <=> $mappings->{$b}} keys(%$mappings);
    my $best_score = $mappings->{$sorted_mappings_by_score[0]};
    $seq_region_2_best_score->{$seq_name} = $best_score;
  }
  # find seq region with highest mapping
  my @seq_regions_with_best_score = ();
  my @sorted_seq_regions_by_score = sort {$seq_region_2_best_score->{$a} <=> $seq_region_2_best_score->{$b}} keys(%$seq_region_2_best_score);
  my $threshold = $seq_region_2_best_score->{$sorted_seq_regions_by_score[0]};
  foreach my $seq_region (keys %$seq_region_2_best_score) {
    if ($seq_region_2_best_score->{$seq_region} <= $threshold) {
      push @seq_regions_with_best_score, $seq_region;
    }
  }
  foreach my $new_chrom (@seq_regions_with_best_score) {
    my $mappings = $filtered->{$new_chrom};
    my @sorted_mappings_by_score = sort {$mappings->{$a} <=> $mappings->{$b}} keys(%$mappings);
    my $best_score = $mappings->{$sorted_mappings_by_score[0]};
    foreach my $mapping (@sorted_mappings_by_score) {
      my $score = $mappings->{$mapping};
      if ($score <= $threshold) {
        $$count_out++;
        my ($new_start, $new_end, $new_strand) = split(':', $mapping);
        push @$output_lines, "$id\t$new_chrom\t$new_start\t$new_end\t$new_strand\t$score\n";
      }
    }
  }
}

sub filter_svf_mapping {
  my $config = shift;
  my $file_prev_mappings = $config->{file_prev_mappings};
  my $file_mappings = $config->{file_mappings};
  my $file_filtered_mappings = $config->{'file_filtered_mappings'};
  my $file_failed_filtered_mappings = $config->{'file_failed_filtered_mappings'};

  my $prev_mappings = _get_prev_mappings($file_prev_mappings);
  my $mapped = _get_svf_mappings($file_mappings, $prev_mappings);

  my $filtered_mappings = {};
  # for each structural variation feature get the best mapping for each coord type
  foreach my $id (keys %$mapped) {
    foreach my $coord_type (keys %{$mapped->{$id}}) {
      my $mappings = $mapped->{$id}->{$coord_type};
      my @sorted_coord_by_score = sort { $mappings->{$b} <=> $mappings->{$a} } keys %$mappings;
      my $algn_score_threshold = 0.5;
      if ($mappings->{$sorted_coord_by_score[0]} == 1.0) {
        $algn_score_threshold = 1.0;
      }
      my $count_exceed_threshold = grep {$_ >= $algn_score_threshold} values %$mappings;
      if ($count_exceed_threshold > 1) {
        my $prev_mapping = $prev_mappings->{$id}->{$coord_type};
        my $mapped_location = best_mapping($algn_score_threshold, $prev_mapping, $mappings);
        $filtered_mappings->{$id}->{$coord_type}->{$mapped_location} = $mappings->{$mapped_location};
      } else { # $count_excess_threshold <= 1 either get 1 mapping or no mapping if threshold to low
        foreach my $coord (@sorted_coord_by_score) {
          my $score = $mappings->{$coord};
          if ($score >= $algn_score_threshold) {
            $filtered_mappings->{$id}->{$coord_type}->{$coord} = $score;
          }
        }
      }
    }
  }
  # check all coord types from last release have been remapped, also check the correct order
  my $final_mappings = {};
  my $fh_failed = FileHandle->new($file_failed_filtered_mappings, 'w');
  my @start_order = qw/outer_start seq_region_start inner_start/;
  my @end_order = qw/inner_end seq_region_end outer_end/;
  foreach my $id (keys %$prev_mappings) {
    my @prev_coord_types = grep {$_ ne 'seq_region_name'} keys %{$prev_mappings->{$id}};
    my @new_coord_types =  keys %{$filtered_mappings->{$id}};
    if (overlap(\@prev_coord_types, \@new_coord_types)) {
      my $start_coords_in_order = coords_are_in_order(\@start_order, $filtered_mappings->{$id});
      my $end_coords_in_order = coords_are_in_order(\@end_order, $filtered_mappings->{$id});
      if (!($start_coords_in_order && $end_coords_in_order)) {
        print $fh_failed "$id Coords not in order\n";
        next;
      }
      # check start coords are smaller than end coords
      my $start = get_start($filtered_mappings->{$id}->{seq_region_start});
      my $end = get_start($filtered_mappings->{$id}->{seq_region_end});
      if ($end < $start) {
        my $swap_map = {
          'outer_start'      => 'inner_end',
          'seq_region_start' => 'seq_region_end',
          'inner_start'      => 'outer_end',

          'inner_end'        => 'outer_start',
          'seq_region_end'   => 'seq_region_start',
          'outer_end'        => 'inner_start',
        };
        my @order = qw/outer_start seq_region_start inner_start inner_end seq_region_end outer_end/;
        my $after_swap_mappings = {};
        foreach my $c (@order) {
          $final_mappings->{$id}->{$c} = get_start($filtered_mappings->{$id}->{$swap_map->{$c}});
        }
        next;
      }
      foreach my $c (@new_coord_types) {
        my $start = get_start($filtered_mappings->{$id}->{$c});
        $final_mappings->{$id}->{$c} = $start;
      }
    } else {
      print $fh_failed "$id Incomplete mappings Prev ", join(', ', @prev_coord_types), " new ", join(', ', @new_coord_types), "\n";
    }
  }
  $fh_failed->close;

  # write filtered mappings to file
  my $fh_filtered_mappings = FileHandle->new($file_filtered_mappings, 'w');
  my $already_stored = {};
  foreach my $id (keys %$final_mappings) {
    my ($svf_id, $strand) = split('_', $id);
    next if ($already_stored->{$svf_id});
    my @values = ();
    my $seq_region_name = $prev_mappings->{$id}->{seq_region_name};
    push @values, "seq_region_name=$seq_region_name";
    push @values, "structural_variation_feature_id=$svf_id";
    push @values, "seq_region_strand=$strand";
    foreach my $coord (keys %{$final_mappings->{$id}}) {
      my $start = $final_mappings->{$id}->{$coord};
      push @values, "$coord=$start";
    }
    print $fh_filtered_mappings join("\t", @values), "\n";
    $already_stored->{$svf_id} = 1;
  }
  $fh_filtered_mappings->close();
}

sub get_start {
  my $coord = shift;
  if ($coord) {
    my @keys = keys %$coord;
    return $keys[0];
  }
  return 0;
}

sub best_mapping {
  my $threshold = shift;
  my $prev_location = shift;
  my $new_mappings = shift;
  my $new_mapping = '';
  my $diff = 1_000_000_000;

  foreach my $start (keys %$new_mappings) {
    my $score = $new_mappings->{$start};
    if ($score >= $threshold) {
      my $new_diff = abs($prev_location - $start);
      if ($diff > $new_diff) {
        $diff = $new_diff;
        $new_mapping = $start;
      }
    }
  }
  return $new_mapping;
}

sub coords_are_in_order {
  my $order = shift;
  my $mappings = shift;

  for (my $i = 0; $i < scalar @$order - 1; $i++) {
    my $a_coord = $order->[$i];
    my $b_coord = $order->[$i + 1];
    my $a_start = get_start($mappings->{$a_coord});
    my $b_start = get_start($mappings->{$b_coord});
    if ($a_start && $b_start) {
      if ($a_start > $b_start) {
        return 0;
      }
    }
  }
  return 1;
}
sub overlap {
  my $a = shift;
  my $b = shift;
  return (scalar @$a == scalar @$b);
}

sub _get_prev_mappings {
  my $file = shift;
  my $prev_mappings = {};
  my $fh = FileHandle->new($file, 'r');
  while (<$fh>) {
    chomp;
    #93742   seq_region_name=X;seq_region_start=44601505;inner_start=44601505;inner_end=44619579;seq_region_end=44619579
    my ($vf_id, $coords) = split/\t/;
    foreach my $coord (split(';', $coords)) {
      my ($key, $value) = split('=', $coord);
      $prev_mappings->{"$vf_id\_1"}->{$key} = $value;
      $prev_mappings->{"$vf_id\_-1"}->{$key} = $value;
    }
  }
  $fh->close();
  return $prev_mappings;
}

sub _get_svf_mappings {
  my $file = shift;
  my $prev_mappings = shift;
  my $fh_mappings = FileHandle->new($file, 'r');
  my $mapped = {};
  while (<$fh_mappings>) {
    chomp;
    #feature_id-coord, coord: outer_start seq_region_start inner_start inner_end seq_region_end outer_end

    #107100:upstream 1       460501  461001  1       6       1       496     0.99001996007984
    #query_name seq_region_name start end strand map_weight edit_dist score_count algn_score
    #102027-seq_region_start 10      44974586        44974758        1       145     1       274     score=0.850746268656716 173M28H
    my ($query_name, $new_seq_name, $new_start, $new_end, $new_strand, $map_weight, $algn_score, $cigar_str) = split("\t", $_);
    my ($id, $coord) = split('-', $query_name);
    $id = "$id\_$new_strand";
    if ($prev_mappings->{$id}->{seq_region_name} eq "$new_seq_name") {
      $mapped->{$id}->{$coord}->{$new_start} = $algn_score;
    }
  }
  $fh_mappings->close();
  return $mapped;
}

sub filter_vf_mapping_first_round {
  my ($config, $stats_config, $chroms) = @_;

  my $use_prior_info = $config->{use_prior_for_filtering};
  my $algn_score_threshold = $config->{algn_score_threshold};
  my $map_to_chrom_only = $config->{map_to_chrom_only};
  my $file_mappings = $config->{file_mappings};
  my $file_filtered_mappings = $config->{file_filtered_mappings};
  my $file_failed_filtered_mappings = $config->{file_failed_filtered_mappings};

  my $multi_map_all = {};
  my $multi_map_same_chrom = {};

  open FH_MAPPINGS, "<$file_mappings" or die $!;
  open FH_FAILED, ">>$file_failed_filtered_mappings" or die $!;
  open FH_FILTERED_MAPPINGS, ">>$file_filtered_mappings" or die $!;

  while (<FH_MAPPINGS>) {
    chomp;
    my ($old_seq_info, $new_seq_info, $query_name, $map_weight, $cigar, $relative_algn_score, $clipped_info) = split("\t", $_);
    #query_name: 156358-150-1-150-11:5502587:5502587:1:T/C:rs202026261:dbSNP:SNV
    #filter old chrom name

    my $old_chrom_name = get_old_chrom_name($query_name) if ($use_prior_info);
    my ($chrom_name, $start, $end, $strand) = split(' ', $new_seq_info);

    if ($relative_algn_score < $algn_score_threshold) {
      $stats_config->{stats_failed_poor_score}++;
      print FH_FAILED "$query_name\t$POOR_SCORE\t$relative_algn_score\n";
    } else {
      if ($map_weight > 1) {
        if ($map_to_chrom_only) {
          $multi_map_all->{$query_name}->{$new_seq_info} = $relative_algn_score if ($chroms->{$chrom_name});
        } else {
          $multi_map_all->{$query_name}->{$new_seq_info} = $relative_algn_score;
        }
        if ($use_prior_info) {
          $multi_map_same_chrom->{$query_name}->{$new_seq_info} = $relative_algn_score if ($chrom_name eq $old_chrom_name);
        }
      } else {
        if ($relative_algn_score >= $algn_score_threshold) {
          my $line = print_feature_line($query_name, $new_seq_info, $relative_algn_score);
          print FH_FILTERED_MAPPINGS $line, "\n";
          $stats_config->{stats_unique_map}++;
        }
      }
    }
  }

  close FH_FAILED;
  close FH_FILTERED_MAPPINGS;
  close FH_MAPPINGS;
  
  return {
    multi_map_all => $multi_map_all,
    multi_map_same_chrom => $multi_map_same_chrom,
    stats_config => $stats_config,
  };

}

sub filter_vf_mapping_second_round {
  my ($config, $stats_config, $multi_map_working) = @_;

  my $algn_score_threshold = $config->{algn_score_threshold};
  my $max_map_weight = $config->{max_map_weight};

  my $file_filtered_mappings = $config->{file_filtered_mappings};
  my $file_failed_filtered_mappings = $config->{file_failed_filtered_mappings};
  open FH_FAILED, ">>$file_failed_filtered_mappings" or die $!;
  open FH_FILTERED_MAPPINGS, ">>$file_filtered_mappings" or die $!;

  my $filtered_multi_map = {};

  foreach my $query_name (keys %$multi_map_working) {
    my $mappings = $multi_map_working->{$query_name};
    my $count_before = scalar keys %$mappings;
    if ($count_before == 0) {
      next; # only maps to non_chrom, can only happen id use_prior_info is set
    }
    my $count_after = 0;

    my @new_seq_infos = sort { $mappings->{$b} <=> $mappings->{$a} } keys(%$mappings);

    # check that the first one has score 1, then only select those with score 1
    # if the first has not score of one then select based on threshold
    my $threshold = $algn_score_threshold;
    if ($mappings->{$new_seq_infos[0]} == 1.0) {
      $threshold = 1.0;
    }
    foreach my $new_seq_info (@new_seq_infos) {
      my $score = $mappings->{$new_seq_info};
      if ($score >= $threshold) {
        my $line = print_feature_line($query_name, $new_seq_info, $score);
        $filtered_multi_map->{$query_name}->{$line} = 1;
        $count_after++;
      }
    }
    if ($count_after == 0) {
      $stats_config->{stats_failed_after_filter}++;
    }
  }

  foreach my $query_name (keys %$filtered_multi_map) {
    my $mappings = $filtered_multi_map->{$query_name};
    my $count = scalar keys %$mappings;
    if ($count >= $max_map_weight) {
      $stats_config->{stats_exceeds_max_map_weight}++;
      print FH_FAILED "$query_name\t$TOO_MANY_MAPPINGS\n";
    } else {
      if ($count == 1) {
        $stats_config->{stats_unique_map_after_filter}++;
      } else {
        $stats_config->{stats_multi_map_after_filter}++;
      }
      foreach my $line (keys %$mappings) {
        print FH_FILTERED_MAPPINGS $line, "\n";
      }
    }
  }
  close FH_FAILED;
  close FH_FILTERED_MAPPINGS;
}

sub print_feature_line {
  my ($query_name, $new_seq_info, $score) = @_;
  my ($seq_name, $start, $end, $strand) = split(' ', $new_seq_info);
  my $line =  join("\t", ($query_name, $seq_name, $start, $end, $strand, $score));
  return $line;
}

sub get_old_chrom_name {
  my $query_name = shift;
  my @query_name_components_fst_part = split('-', $query_name, 5);
  my @query_name_components_snd_part = split(':', $query_name_components_fst_part[4], 2);
  return $query_name_components_snd_part[0];
}

sub qc_mapped_vf {
  my $config = shift;
  my $mapped_features_file = $config->{mapped_features_file};
  my $update_features_file = $config->{update_features_file};
  my $failure_reasons_file = $config->{failure_reasons_file};  
  my $feature_table = $config->{feature_table};
  my $vdba = $config->{vdba};
  my $vdba_oldasm = $config->{vdba_oldasm};
  my $fasta_db = Bio::DB::Fasta->new($config->{fasta_db});
  my @ids      = $fasta_db->get_all_primary_ids;
  my $sequence_name_2_id = {};
  foreach my $id (@ids) {
    if ($id =~ /:/) {
      my @components = split/:/, $id;
      my $sequence_name = $components[2];
      $sequence_name_2_id->{$sequence_name} = $id;
    } else {
      $sequence_name_2_id->{$id} = $id;
    }
  }

  my $fh = FileHandle->new($mapped_features_file, 'r');
  my $fh_update = FileHandle->new($update_features_file, 'w');
  my $fh_failure_reasons = FileHandle->new($failure_reasons_file, 'w');

  my $failed_variants_oldasm = failed_variations($vdba_oldasm);
  my $failed_variants_newasm = {};

  while (<$fh>)  {
    chomp;
    my $data = read_line($_);
    my $vf_id = $data->{variation_feature_id};
    my $variation_id = $data->{variation_id};
    my $variation_name = $data->{variation_name};
    my $seq_region_id = $data->{seq_region_id};
    my $start = $data->{seq_region_start};
    my $end = $data->{seq_region_end};
    my $seq_region_name = $data->{seq_region_name};
    my $seq_region_strand = $data->{seq_region_strand};
    my $map_weight = $data->{map_weight};
    my $allele_string = $data->{allele_string};
    my $alignment_quality = $data->{alignment_quality};

    my @old_failure_reasons = keys %{$failed_variants_oldasm->{$variation_id}};
 # To check after remapping variants to a new assembly:
#    1 =>  'Variant maps to more than 3 different locations',
#    2 =>  'None of the variant alleles match the reference allele',
#    5 =>  'Variant does not map to the genome',
#    15 => 'Mapped position is not compatible with reported alleles',
#    17 => 'Variant can not be re-mapped to the current assembly',
#    19 => 'Variant maps to more than one genomic location',

    foreach my $failed_id (@old_failure_reasons) {
      $failed_variants_newasm->{$variation_id}->{$failed_id} = 1;
    }
    if ($map_weight > 1) {
      $failed_variants_newasm->{$variation_id}->{19} = 1;
    }
    if ($seq_region_strand == -1) {
      print $fh_update "UPDATE $feature_table SET seq_region_strand=1 WHERE variation_feature_id=$vf_id;\n";
    }

    if ($end >= $start) {
      my @allele_string_init = sort split('/', $allele_string);
      my @allele_string_rev_comp = sort split ('/', $allele_string);
      foreach my $allele (@allele_string_rev_comp) {
        reverse_comp(\$allele);
      }
      my $id = $sequence_name_2_id->{$seq_region_name};
      my $ref =  uc $fasta_db->seq("$id:$start,$end");
      my @new_allele_strings = ();
      my $is_reverse_comp = 0;
      if ( grep( /^$ref$/, @allele_string_init ) ) {
        push @new_allele_strings, $ref;
        foreach my $allele (@allele_string_init) {
          if ($allele ne $ref) {
            push @new_allele_strings, $allele
          }
        }
      } elsif ( grep( /^$ref$/, @allele_string_rev_comp) )  {
        $is_reverse_comp = 1;
        push @new_allele_strings, $ref;
        foreach my $allele (@allele_string_rev_comp) {
          if ($allele ne $ref) {
            push @new_allele_strings, $allele
          }
        }
      } else {
        # 2 None of the variant alleles match the reference allele
        $failed_variants_newasm->{$variation_id}->{2} = 1;
      }
      if (scalar @new_allele_strings > 0) {
        my $new_allele_string = join('/', @new_allele_strings);
        if ($allele_string ne $new_allele_string) {
          print $fh_update "UPDATE $feature_table SET allele_string='$new_allele_string' WHERE variation_feature_id=$vf_id;\n";
        }
        if ($is_reverse_comp) {
          print $fh_update "UPDATE $feature_table SET flip=1 WHERE variation_feature_id=$vf_id;\n";
        }
      }
    } else {
      if ($allele_string =~ /-/) {
        # resort order of alleles with dash first
        my @new_allele_strings = ('-');
        foreach my $allele (split/\//, $allele_string) {
          if ($allele ne '-') {
            push @new_allele_strings, $allele;
          }
        }
        my $new_allele_string = join('/', @new_allele_strings);
        print $fh_update "UPDATE $feature_table SET allele_string='$new_allele_string' WHERE variation_feature_id=$vf_id;\n";
      } else {
        # Mapped position is not compatible with reported alleles 15
        $failed_variants_newasm->{$variation_id}->{15} = 1;
      }
    }
  }

  # write failure reasons

  $fh->close();
  $fh_update->close();
  foreach my $variant_id (keys %$failed_variants_newasm) {
    print $fh_failure_reasons $variant_id, "\t", join("\t", sort keys %{$failed_variants_newasm->{$variant_id}}), "\n";
  }
  $fh_failure_reasons->close();
}

sub read_line {
  my $line = shift;
  my @key_values = split("\t", $line);
  my $mapping = {};
  foreach my $key_value (@key_values) {
    my ($table_name, $value) = split('=', $key_value, 2);
    $mapping->{$table_name} = $value;
  }
  return $mapping;
}

sub failed_variations {
  my $vdba = shift;
  my $dbh = $vdba->dbc->db_handle;
  my $failed_variations = {};
  my $sth = $dbh->prepare(qq{
    SELECT variation_id, failed_description_id FROM failed_variation WHERE failed_description_id not in (1, 2, 15, 17, 19);
  }, {mysql_use_result => 1});
  $sth->execute();
  while (my $row = $sth->fetchrow_arrayref) {
    $failed_variations->{$row->[0]}->{$row->[1]} = 1;
  }
  $sth->finish();
  return $failed_variations;
}


1;
