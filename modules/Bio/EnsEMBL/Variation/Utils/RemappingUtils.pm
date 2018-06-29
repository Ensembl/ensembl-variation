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

our @EXPORT_OK = qw(map_variant filter_vf_mapping_first_round filter_vf_mapping_second_round qc_mapped_vf);

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
        print $fh_update "UPDATE $feature_table SET allele_string='$new_allele_string' WHERE variation_feature_id=$vf_id;\n";
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
