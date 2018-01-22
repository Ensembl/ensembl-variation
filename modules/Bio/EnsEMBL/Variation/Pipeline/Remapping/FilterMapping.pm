#!/usr/bin/env perl
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




=head1 CONTACT

Please email comments or questions to the public Ensembl
developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

Questions may also be sent to the Ensembl help desk at
<helpdesk.org>.

=cut
package Bio::EnsEMBL::Variation::Pipeline::Remapping::FilterMapping;

use strict;
use warnings;

use FileHandle;
use Bio::DB::Fasta;
use Bio::EnsEMBL::Registry;

use base ('Bio::EnsEMBL::Hive::Process');


# failure_reasons
my $NO_MAPPING = 'no mapping';
my $TOO_MANY_MAPPINGS = 'map weight > 5';
my $POOR_SCORE = 'does not pass alignment score';

sub fetch_input {
  my $self = shift;

  # initialise file names
  my $params = {
    mapping_results_dir   => 'file_mappings',
    filtered_mappings_dir => 'file_filtered_mappings',
    statistics_dir        => 'file_statistics',
    dump_features_dir     => 'file_init_feature',
    load_features_dir     => 'file_load_features',
    fasta_files_dir       => 'fasta_file',
  };
  my $file_number = $self->param('file_number');
  my $individual_dir = '';
  if ($self->param('mode') eq 'remap_read_coverage') {
    my $individual_id = $self->param('individual_id');
    $individual_dir = "/$individual_id/";
  }
  foreach my $param (keys %$params) {
    my $dir = $self->param($param) . $individual_dir;
    if ($param =~ /mapping_results_dir/) {
      my $file_mappings = "$dir/mappings_$file_number.txt";
      my $file_failed_mappings = "$dir/failed_mapping_$file_number.txt";
      $self->param('file_mappings', $file_mappings);
      $self->param('file_failed_mappings', $file_failed_mappings);
    } elsif ($param =~ /filtered_mappings_dir/) {
      my $file_filtered_mappings = "$dir/$file_number.txt";
      my $file_failed_filtered_mappings = "$dir/failed_mappings_$file_number.txt";
      $self->param('file_filtered_mappings', $file_filtered_mappings);
      $self->param('file_failed_filtered_mappings', $file_failed_filtered_mappings);
    } elsif ($param =~ /fasta_files_dir/) {
      my $fasta_file = "$dir/$file_number.fa";
      $self->param('fasta_file', $fasta_file);
    } else {
      my $file = "$dir/$file_number.txt";
      my $file_param = $params->{$param};
      $self->param($file_param, $file);
    }
  }

  # get seq_region ids for new assembly
  my $registry = 'Bio::EnsEMBL::Registry';
  $registry->load_all($self->param('registry_file_newasm'));
  my $cdba = $registry->get_DBAdaptor($self->param('species'), 'core');
  $self->param('cdba', $cdba);
  my $vdba = $registry->get_DBAdaptor($self->param('species'), 'variation');
  $self->param('vdba', $vdba);

}

sub run {
  my $self = shift;
  my $mode = $self->param('mode');
  if ($mode eq 'remap_multi_map') {
    $self->report_failed_mappings();
    $self->filter_mapping_results_dbsnp();
    $self->join_feature_data();
  } elsif ($mode eq 'remap_alt_loci') {
    $self->filter_mapping_results_alt_loci();
    $self->report_failed_mappings();
    $self->join_feature_data();
  } elsif ($mode eq 'remap_read_coverage') {
    $self->report_failed_read_coverage_mappings();
    $self->filter_read_coverage_mapping_results();
    $self->join_read_coverage_data();
  } elsif ($mode eq 'remap_svf') {
    $self->filter_svf_mapping_results();
    $self->join_svf_data();
  } elsif ($mode eq 'remap_qtls') {
    $self->report_failed_read_coverage_mappings();
    $self->filter_read_coverage_mapping_results();
    $self->join_qtl_data();
  } else {
    $self->report_failed_mappings();
    $self->filter_mapping_results();
    $self->join_feature_data();
  }
  $self->write_statistics();

}

sub write_output {
  my $self = shift;
}

sub report_failed_read_coverage_mappings {
  my $self = shift;
  my $count_mapped = 0;
  my $count_unmapped = 0;
  my $file_mappings = $self->param('file_mappings');
  my $file_failed_mappings = $self->param('file_failed_mappings');

  my $mapped = {};
  my $unmapped = {};
  my $fh_mappings = FileHandle->new($file_mappings, 'r');
  while (<$fh_mappings>) {
    chomp;
    my ($query_name, $seq_name, $start, $end, $strand, $map_weight, $score) = split("\t", $_);
    my ($id, $type) = split(':', $query_name);
    $mapped->{$id} = 1;
  }
  $fh_mappings->close();

  $count_mapped = scalar keys %$mapped;

  $self->param('pre_count_mapped', $count_mapped);
  $self->param('pre_count_unmapped', $count_unmapped);
  my $count_input_ids = 0;
  my $fasta_file = $self->param('fasta_file');
  my $fh_fasta_file = FileHandle->new($fasta_file, 'r');
  while (<$fh_fasta_file>) {
    chomp;
    if (/^>/) {
      $count_input_ids++;
    }
  }
  $fh_fasta_file->close();
  $self->param('count_input_ids', $count_input_ids);
}

sub report_failed_mappings {
  my $self = shift;
  my $file_mappings = $self->param('file_mappings');
  my $file_failed_mappings = $self->param('file_failed_mappings');

  my $mapped = {};
  my $unmapped = {};
  my $fh_mappings = FileHandle->new($file_mappings, 'r');
  while (<$fh_mappings>) {
    chomp;
    my ($old_seq_info, $new_seq_info, $query_name, $map_weight, $cigar, $relative_algn_score, $clipped_info) = split("\t", $_);
    $mapped->{$query_name} = 1;
  }
  $fh_mappings->close();
  
  my $fh_failed_mappings = FileHandle->new($file_failed_mappings, 'r');
  while (<$fh_failed_mappings>) {
    chomp;
    my ($indel, $old_seq_info, $new_seq_info, $query_name, $map_weight, $cigar, $relative_algn_score, $clipped_info) = split("\t", $_);
    unless ($mapped->{$query_name}) {
      $unmapped->{$query_name} = 1;
    }
  }
  $fh_failed_mappings->close();

  my $file_failed_filtered_mappings = $self->param('file_failed_filtered_mappings');
  my $fh_failed_filtered_mappings = FileHandle->new($file_failed_filtered_mappings, 'w');
  foreach my $query_name (keys %$unmapped) {
    print $fh_failed_filtered_mappings "$query_name\t$NO_MAPPING\n";
  }
  $fh_failed_filtered_mappings->close();

  my $count_mapped = scalar keys %$mapped;
  my $count_unmapped = scalar keys %$unmapped;

  $self->param('pre_count_mapped', $count_mapped);
  $self->param('pre_count_unmapped', $count_unmapped);

  # get input id counts
  my $count_input_ids = 0;
  my $fasta_file = $self->param('fasta_file');
  my $fh_fasta_file = FileHandle->new($fasta_file, 'r');
  while (<$fh_fasta_file>) {
    chomp;
    if (/^>/) {
      $count_input_ids++;
    }
  }
  $fh_fasta_file->close();
  $self->param('count_input_ids', $count_input_ids);

}

sub get_seq_region_names {
  my $self = shift;
  my $coord = shift;
  my $cdba = $self->param('cdba');  
  my $seq_region_names = {};
  my $sa = $cdba->get_SliceAdaptor;
  foreach my $slice (@{$sa->fetch_all($coord)}) {
    $seq_region_names->{$slice->seq_region_name} = 1;
  } 
  return $seq_region_names; 
}

sub get_old_chrom_name {
  my $self = shift;
  my $query_name = shift;
  my @query_name_components_fst_part = split('-', $query_name, 5);
  my @query_name_components_snd_part = split(':', $query_name_components_fst_part[4], 2);
  return $query_name_components_snd_part[0];
}

sub filter_mapping_results {
  my $self = shift;

  my $algn_score_threshold = $self->param('algn_score_threshold');
  my $use_prior_info = $self->param('use_prior_for_filtering');
  my $map_to_chrom_only = $self->param('map_to_chrom_only');
  my $max_map_weight = $self->param('max_map_weight');
  my $chroms = $self->get_seq_region_names('chromosome') if ($map_to_chrom_only);

  my $file_mappings = $self->param('file_mappings');
  my $fh_mappings = FileHandle->new($file_mappings, 'r');

  my $file_filtered_mappings = $self->param('file_filtered_mappings');
  my $fh_filtered_mappings = FileHandle->new($file_filtered_mappings, 'w');

  my $file_failed_filtered_mappings = $self->param('file_failed_filtered_mappings');

  open FH_FAILED, ">>$file_failed_filtered_mappings" or die $!;

  my $multi_map_all = {};
  my $multi_map_same_chrom = {};

  my ($stats_unique_map, $stats_multi_map, $stats_failed_poor_score, $stats_unique_map_after_filter, $stats_multi_map_after_filter, $stats_failed_after_filter, $stats_failed_non_chrom, $stats_exceeds_max_map_weight);

  while (<$fh_mappings>) {
    chomp;
    my ($old_seq_info, $new_seq_info, $query_name, $map_weight, $cigar, $relative_algn_score, $clipped_info) = split("\t", $_);
    #query_name: 156358-150-1-150-11:5502587:5502587:1:T/C:rs202026261:dbSNP:SNV
    #filter old chrom name

    my $old_chrom_name = $self->get_old_chrom_name($query_name) if ($use_prior_info);    

    my ($chrom_name, $start, $end, $strand) = split(' ', $new_seq_info);
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
        my $line = $self->print_feature_line($query_name, $new_seq_info, $relative_algn_score);
        print $fh_filtered_mappings $line, "\n";
        $stats_unique_map++;
      } else {
        $stats_failed_poor_score++;
        print FH_FAILED "$query_name\t$POOR_SCORE\n";
      }
    }
  }
  $fh_mappings->close();

  my $count_multi_map_all        = scalar keys %$multi_map_all;
  my $count_multi_map_same_chrom = scalar keys %$multi_map_same_chrom;

  $stats_failed_non_chrom = $count_multi_map_all - $count_multi_map_same_chrom;

  my $multi_map_working = ($use_prior_info) ? $multi_map_same_chrom : $multi_map_all;
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
        my $line = $self->print_feature_line($query_name, $new_seq_info, $score);
        $filtered_multi_map->{$query_name}->{$line} = 1;
        $count_after++;
      }
    }
    if ($count_after == 0) {
      $stats_failed_after_filter++;
    }
  }

  foreach my $query_name (keys %$filtered_multi_map) {
    my $mappings = $filtered_multi_map->{$query_name};
    my $count = scalar keys %$mappings;
    if ($count >= $max_map_weight) {
      $stats_exceeds_max_map_weight++;
      print FH_FAILED "$query_name\t$TOO_MANY_MAPPINGS\n";
    } else {
      if ($count == 1) {
        $stats_unique_map_after_filter++;
      } else {
        $stats_multi_map_after_filter++;
      }
      foreach my $line (keys %$mappings) {
        print $fh_filtered_mappings $line, "\n";
      }
    }
  }
  $self->param('stats_unique_map', $stats_unique_map);
  $self->param('stats_unique_map_after_filter', $stats_unique_map_after_filter);
  $self->param('stats_multi_map', $stats_multi_map);
  $self->param('stats_multi_map_after_filter', $stats_multi_map_after_filter);
  $self->param('stats_failed_poor_score', $stats_failed_poor_score);
  $self->param('stats_failed_after_filter', $stats_failed_after_filter);
  $self->param('stats_failed_non_chrom', $stats_failed_non_chrom);
  $self->param('stats_exceeds_max_map_weight', $stats_exceeds_max_map_weight);

  $fh_filtered_mappings->close();
  close FH_FAILED;
}

sub filter_mapping_results_dbsnp {
  my $self = shift;
  my $algn_score_threshold = $self->param('algn_score_threshold');

  my $file_mappings = $self->param('file_mappings');
  my $fh_mappings = FileHandle->new($file_mappings, 'r');

  my $file_filtered_mappings = $self->param('file_filtered_mappings');
  my $fh_filtered_mappings = FileHandle->new($file_filtered_mappings, 'w');

  my ($stats_failed, $stats_unique_map, $stats_multi_map);

  my $mapped = {};
  while (<$fh_mappings>) {
    chomp;
    # old_seq_info is not in available
    #       1 19432 19432 1 480394.0-101-1-101-A/G:rs75513536       8       203M    1                       clipped nucleotides 0
    #       1 19432 19432 1 480394.1-101-1-101-A/G:rs75513536       8       203M    0.995073891625616       clipped nucleotides 0

    #        1 18849 18849 -1        14086.1-110-1-497-G:C/G:rs708635        3       41S567M 0.932565789473684       clipped nucleotides 41
    #        1 19172 19172 -1        14567.1-200-1-200-T:C/T:rs806720        8       110M2I289M      0.992518703241895       clipped nucleotides 0
    #        1 19241 19241 -1        4257.1-257-1-250-G:A/G:rs380444 3       508M    0.998031496062992       clipped nucleotides 0
    my ($old_seq_info, $new_seq_info, $query_name, $map_weight, $cigar, $relative_algn_score, $clipped_info) = split("\t", $_);

    my ($chrom_name, $start, $end, $strand) = split(' ', $new_seq_info);
    my @query_name_components = split('-', $query_name, 5);
    my $vf_id_info = $query_name_components[0];
    my ($vf_id, $version) = split('\.', $vf_id_info);
    my ($allele, $allele_string, $rsid) = split(':', $query_name_components[4]);
    $mapped->{$vf_id}->{$new_seq_info}->{$query_name} = $relative_algn_score;
  }
  $fh_mappings->close();

  #content of $mapped could look like this:
  #3509 22 18401349 18401376 -1 3509.1-200-28-200--/CGGAGCCAGAGGGCCGGGGGGTCCCACA:rs361903 0.992990654205608
  #3509 22 18718507 18718506 1 3509.0-200-0-200--/CGGAGCCAGAGGGCCGGGGGGTCCCACA:rs361903 1
  #3509 22 18223193 18223220 1 3509.1-200-28-200--/CGGAGCCAGAGGGCCGGGGGGTCCCACA:rs361903 0.992990654205608
  #3509 22 21344289 21344316 1 3509.1-200-28-200--/CGGAGCCAGAGGGCCGGGGGGTCCCACA:rs361903 0.995327102803738
  #3509 22 18650965 18650992 -1 3509.1-200-28-200--/CGGAGCCAGAGGGCCGGGGGGTCCCACA:rs361903 0.995327102803738

  #27587 9 65226786 65226786 -1 27587.1-140-1-507-C/G:rs1838599 0.998456790123457
  #27587 9 65226786 65226786 -1 27587.0-140-1-507-C/G:rs1838599 1

  #27587 9 63022076 63022076 1 27587.1-140-1-507-C/G:rs1838599 0.970679012345679
  #27587 9 63022076 63022076 1 27587.0-140-1-507-C/G:rs1838599 0.972222222222222

  #27587 9 65795203 65795203 1 27587.1-140-1-507-C/G:rs1838599 0.976851851851852
  #27587 9 65795203 65795203 1 27587.0-140-1-507-C/G:rs1838599 0.978395061728395

  #27587 9 64741286 64741286 -1 27587.1-140-1-507-C/G:rs1838599 0.983024691358025
  #27587 9 64741286 64741286 -1 27587.0-140-1-507-C/G:rs1838599 0.984567901234568

  # if coords are the same, choose the one with better score --> consider allele string?

  foreach my $vf_id (keys %$mapped) {
    my $passed = {};

    foreach my $new_seq_info (keys %{$mapped->{$vf_id}}) {
    # query_name to score
      my $query_name_to_score = $mapped->{$vf_id}->{$new_seq_info};
      my @scores = sort { $query_name_to_score->{$b} <=> $query_name_to_score->{$a} } keys(%$query_name_to_score);
      my $query_name = $scores[0];
      my $score = $query_name_to_score->{$query_name};

      if ($score >= $algn_score_threshold) {
        $passed->{"$vf_id\t$new_seq_info\t$query_name"} = $score;
      }
    }

    my $count_passed = scalar keys %$passed;
    if ($count_passed > 0) {
    # filter by score
      my $filtered = {};
      my @keys = sort {$passed->{$b} <=> $passed->{$a}} keys(%$passed);
      my $threshold = $algn_score_threshold;
      if ($passed->{$keys[0]} == 1.0) {
        $threshold = 1.0;
      }
      foreach my $result (keys %$passed) {
        my $score = $passed->{$result};
        if ($score >= $threshold) {
          $filtered->{$result} = $score;
          my ($vf_id, $new_seq_info, $query_name) = split("\t", $result);
          # remove version for vf_id: 
          $query_name =~ s/\.\d+//g;
          my $line = $self->print_feature_line($query_name, $new_seq_info, $score);
          print $fh_filtered_mappings $line, "\n";
        }   
      }

      $count_passed = scalar keys %$filtered;
      if ($count_passed == 0) {
        $stats_failed++;
      } elsif ($count_passed == 1) {
        $stats_unique_map++;
      } else {
        $stats_multi_map++;
      }
    } else {
      $stats_failed++;
    }
  }

  $fh_filtered_mappings->close();

  $self->param('stats_unique_map', $stats_unique_map);
  $self->param('stats_multi_map', $stats_multi_map);
  $self->param('stats_failed', $stats_failed);
}

sub filter_mapping_results_alt_loci {
  my $self = shift;

  my $alt_loci_to_ref = {};
  my $alt_loci_coords = {};
  my $ref_to_unique_region_coords = {}; 
  my $cdba = $self->param('cdba');
  my $sa = $cdba->get_SliceAdaptor;
  my $aefa = $cdba->get_AssemblyExceptionFeatureAdaptor;
  my $slices = $sa->fetch_all('chromosome', undef, 1, 1);

  foreach my $slice (@$slices) {
    my $seq_region_name = $slice->seq_region_name;
    my $assembly_exception_type = $slice->assembly_exception_type; # 'HAP', 'PAR', 'PATCH_FIX' or 'PATCH_NOVEL'
      if ($assembly_exception_type eq 'HAP') {
        my $alt_slices = $sa->fetch_by_region_unique('chromosome', $seq_region_name);
        foreach my $alt_slice (@$alt_slices) {
          my $start = $alt_slice->start;
          my $end = $alt_slice->end;
          $alt_loci_coords->{$seq_region_name}->{start} = $start;
          $alt_loci_coords->{$seq_region_name}->{end} = $end;
        }
      } elsif ($assembly_exception_type eq 'REF') {
        my $assembly_exception_features = $aefa->fetch_all_by_Slice($slice);
        my $ref_start = $slice->start;
        my $ref_end = $slice->end;
        $ref_to_unique_region_coords->{start} = $ref_end;
        $ref_to_unique_region_coords->{end} = $ref_start;
        foreach my $feature (@$assembly_exception_features) {
          my $alt_slice = $feature->alternate_slice();
          my $alt_slice_name = $alt_slice->seq_region_name;
          unless ($alt_slice_name =~ /^X$|^Y$|PATCH/) {
            $alt_loci_to_ref->{$alt_slice_name}->{$seq_region_name} = 1;
          }
        }
      }
  }

  my $algn_score_threshold = $self->param('algn_score_threshold');

  my $file_mappings = $self->param('file_mappings');
  my $fh_mappings = FileHandle->new($file_mappings, 'r');

  my $file_filtered_mappings = $self->param('file_filtered_mappings');
  my $fh_filtered_mappings = FileHandle->new($file_filtered_mappings, 'w');

  my ($stats_failed, $stats_unique_map, $stats_multi_map);

  my $mappings = {};
  while (<$fh_mappings>) {
    chomp;
    my ($old_seq_info, $new_seq_info, $query_name, $map_weight, $cigar, $relative_algn_score, $clipped_info) = split("\t", $_);
    #        CHR_HSCHR1_1_CTG3 118061 118060 1       160136-150-0-150-11:5620652:5620651:1:-/ATTT:rs72401051:dbSNP:insertion 107     112H165M23H     0.503333333333333       clipped nucleotides 135
    #        CHR_HSCHR1_1_CTG3 118082 118082 1       160137-150-1-150-11:5620673:5620673:1:G/A:rs375926866:dbSNP:SNV 112     91H165M45H      0.501661129568106       clipped nucleotides 136
    #        CHR_HSCHR1_1_CTG3 118144 118144 1       160139-150-1-150-11:5620735:5620735:1:G/A:rs7952456:dbSNP:SNV   126     29H165M107H     0.501661129568106       clipped nucleotides 136
    #        CHR_HSCHR1_1_CTG3 118182 118182 1       175371-150-1-150-11:6171720:6171720:1:A/G:rs11040769:dbSNP:SNV  189     11H167M1D66M2I53M2H     0.81063122923588        clipped nucleotides 13

    # query_name: 156358-150-1-150-11:5502587:5502587:1:T/C:rs202026261:dbSNP:SNV

    # filter old chrom name
    my @query_name_components_fst_part = split('-', $query_name, 5);
    my @query_name_components_snd_part = split(':', $query_name_components_fst_part[4], 2);
    my $old_chrom_name = $query_name_components_snd_part[0];

    my ($alt_loci_name, $start, $end, $strand) = split(' ', $new_seq_info);
    if ($alt_loci_to_ref->{$alt_loci_name}->{$old_chrom_name}) {
      my $alt_loci_start = $alt_loci_coords->{$alt_loci_name}->{start};
      my $alt_loci_end   = $alt_loci_coords->{$alt_loci_name}->{end};
      my $updated_start  = $alt_loci_start + $start - 1;
      my $updated_end    = $alt_loci_start + $end - 1;
      my $updated_seq_info = join(' ', $alt_loci_name, $updated_start, $updated_end, $strand);
      # check mappings:
      if ($start < 150 || ($end > ($alt_loci_end - 150) ) ) {
        $mappings->{$query_name}->{$alt_loci_name}->{$updated_seq_info} = $relative_algn_score;
      #print STDERR "Edge mapping $_\n";
      } else {
        if ($relative_algn_score >= $algn_score_threshold ) {
          $mappings->{$query_name}->{$alt_loci_name}->{$updated_seq_info} = $relative_algn_score;
        }
      }
    }
  }
  $fh_mappings->close();

  foreach my $query_name (keys %$mappings) {
    foreach my $chrom (keys %{$mappings->{$query_name}}) {
      foreach my $new_seq_info (keys %{$mappings->{$query_name}->{$chrom}}) {
        my $score = $mappings->{$query_name}->{$chrom}->{$new_seq_info};
        my $line = $self->print_feature_line($query_name, $new_seq_info, $score);
        print $fh_filtered_mappings $line, "\n";
      }
    }
  }
  $fh_filtered_mappings->close();

  $self->param('stats_unique_map', $stats_unique_map);
  $self->param('stats_multi_map', $stats_multi_map);
  $self->param('stats_failed', $stats_failed);
}

sub filter_svf_mapping_results {
  my $self = shift;
  my $algn_score_threshold = $self->param('algn_score_threshold');

  my $file_mappings = $self->param('file_mappings');
  my $fh_mappings = FileHandle->new($file_mappings, 'r');

  my $lookup = {};
  my $file_number = $self->param('file_number');
  my $dump_features_dir = $self->param('dump_features_dir'); 
  my $fh_lookup = FileHandle->new("$dump_features_dir/lookup_$file_number.txt", 'r');
  while (<$fh_lookup>) {
    chomp;
    #93742   seq_region_name=X;seq_region_start=44601505;inner_start=44601505;inner_end=44619579;seq_region_end=44619579 
    my ($vf_id, $coords) = split/\t/;
    foreach my $coord (split(';', $coords)) {
      my ($key, $value) = split('=', $coord);
      $lookup->{$vf_id}->{$key} = $value;
    }
  }
  $fh_lookup->close();
  

  my ($stats_failed, $stats_unique_map, $stats_multi_map);

  my $mapped = {};
  while (<$fh_mappings>) {
    chomp;

    #feature_id-coord, coord: outer_start seq_region_start inner_start inner_end seq_region_end outer_end

    #107100:upstream 1       460501  461001  1       6       1       496     0.99001996007984   
    #query_name seq_region_name start end strand map_weight edit_dist score_count algn_score
    #102027-seq_region_start 10      44974586        44974758        1       145     1       274     score=0.850746268656716 173M28H
    my ($query_name, $new_seq_name, $new_start, $new_end, $new_strand, $map_weight, $algn_score, $cigar_str) = split("\t", $_); 

    my ($id, $coord) = split('-', $query_name);   
    if ($lookup->{$id}->{seq_region_name} eq "$new_seq_name") {
      $mapped->{$id}->{$coord}->{$new_start} = $algn_score;
    }
  }
  $fh_mappings->close();

  my $filtered_mappings = {};
  foreach my $id (keys %$mapped) {
    foreach my $coord_type (keys %{$mapped->{$id}}) {
      my $mappings = $mapped->{$id}->{$coord_type};
      my @coord_2_scores = sort { $mappings->{$b} <=> $mappings->{$a} } keys %$mappings;
      my $algn_score_threshold = 0.75;
      if ($mappings->{$coord_2_scores[0]} == 1.0) {
        $algn_score_threshold = 1.0;
      }
      my $count_exceed_threshold = grep {$_ >= $algn_score_threshold} values %$mappings;

      if ($count_exceed_threshold > 1) {
        my $prev_mapping = $lookup->{$id}->{$coord_type};
        my $mapped_location = best_mapping($algn_score_threshold, $prev_mapping, $mappings);
        $filtered_mappings->{$id}->{$coord_type}->{$mapped_location} = $mappings->{$mapped_location};
      } else {
        foreach my $coord (@coord_2_scores) {
          my $score = $mappings->{$coord};
          if ($score >= $algn_score_threshold) {
            $filtered_mappings->{$id}->{$coord_type}->{$coord} = $score;
          }
        }
      }
    }
  }

  my $final_mappings = {};
  my $failed_mappings = {};
  foreach my $id (keys %$lookup) {
    my @prev_coord_types = grep {$_ ne 'seq_region_name'} keys %{$lookup->{$id}};
    my @new_coord_types =  keys %{$filtered_mappings->{$id}};
    if (overlap(\@prev_coord_types, \@new_coord_types)) {
      my @start_order = qw/outer_start seq_region_start inner_start/;
      my @end_order = qw/inner_end seq_region_end outer_end/;
      my $start_coords_in_order = coords_are_in_order(\@start_order, $filtered_mappings->{$id});
      my $end_coords_in_order = coords_are_in_order(\@end_order, $filtered_mappings->{$id});
      if (!($start_coords_in_order && $end_coords_in_order)) {
        $failed_mappings->{$id} = 'Coords not in order';
        next;
      }
      # check start coords are smaller that end coords
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
      $failed_mappings->{$id} = 'Incomplete mappings';
    }
  }

  my $file_filtered_mappings = $self->param('file_filtered_mappings');
  my $fh_filtered_mappings = FileHandle->new($file_filtered_mappings, 'w');

  foreach my $id (keys %$final_mappings) {
    my @values = ();
    my $seq_region_name = $lookup->{$id}->{seq_region_name};
    push @values, "seq_region_name=$seq_region_name";
    push @values, "structural_variation_feature_id=$id";
    foreach my $coord (keys %{$final_mappings->{$id}}) {
      my $start = $final_mappings->{$id}->{$coord};
      push @values, "$coord=$start";
    }
    print $fh_filtered_mappings join("\t", @values), "\n";
  }

  $fh_filtered_mappings->close();
}

sub filter_read_coverage_mapping_results {
  my $self = shift;
  my $algn_score_threshold = $self->param('algn_score_threshold');

  my $file_mappings = $self->param('file_mappings');
  my $fh_mappings = FileHandle->new($file_mappings, 'r');

  my ($stats_failed, $stats_unique_map, $stats_multi_map);

  #my $fh = FileHandle->new('/lustre/scratch110/ensembl/at7/release_86/chicken/remapping_qtls/debug_filter_mappings', 'w');

  my $mapped = {};
  while (<$fh_mappings>) {
    chomp;
    #107100:upstream 1       460501  461001  1       6       1       496     0.99001996007984   
    #query_name seq_region_name start end strand map_weight edit_dist score_count algn_score
    #5844:upstream   5       30870168        30870668        1       1       1       501M
    my ($query_name, $new_seq_name, $new_start, $new_end, $new_strand, $map_weight, $algn_score, $cigar) = split("\t", $_); 

    my ($id, $type) = split(':', $query_name);   

    $mapped->{$id}->{$type}->{$new_seq_name}->{"$new_start:$new_end:$new_strand"} = $algn_score;
   # print $fh "$id $type $new_seq_name $new_start $new_end $new_strand $algn_score\n";
  }
#  $fh->close;
  $fh_mappings->close();

  my $file_init_feature = $self->param('file_init_feature');
  my $fh_init_feature = FileHandle->new($file_init_feature, 'r'); 

  my $feature_data = {};
  while (<$fh_init_feature>) {
    chomp;
    my $data = $self->read_data($_);
    $feature_data->{$data->{entry}} = $data;
  }
  $fh_init_feature->close();

  my $file_filtered_mappings = $self->param('file_filtered_mappings');
  my $fh_filtered_mappings = FileHandle->new($file_filtered_mappings, 'w');

  foreach my $id (keys %$mapped) {
    my $init_data = $feature_data->{$id};
    my $old_chrom = $init_data->{seq_region_name};
    my $old_start = $init_data->{seq_region_start}; 
    my $old_end = $init_data->{seq_region_end};
    my $diff = $old_end - $old_start + 1;
    my $types = scalar keys %{$mapped->{$id}};
    my $count_out = 0;
    my @output_lines = ();
    if ($types == 2) {
      my $filtered = {};
      foreach my $seq_name (keys %{$mapped->{$id}->{upstream}}) {
        my $filtered_upstream   = $mapped->{$id}->{upstream}->{$seq_name};
        my $filtered_downstream = $mapped->{$id}->{downstream}->{$seq_name};
        next unless ($filtered_upstream && $filtered_downstream);
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
#            $filtered->{$seq_name}->{"$up_start:$down_end:$down_strand"} = $diff_score;
          }
        }
      }           
      if ($filtered->{$old_chrom}) {
        my $mappings = $filtered->{$old_chrom};
        my @keys = sort {$mappings->{$a} <=> $mappings->{$b}} keys(%$mappings);
        my $threshold = $mappings->{$keys[0]};
        foreach my $mapping (@keys) {
          my $score = $mappings->{$mapping};
          if ($score <= $threshold) {
            $count_out++;
            my ($new_start, $new_end, $new_strand) = split(':', $mapping);
            
            push @output_lines, "$id\t$old_chrom\t$new_start\t$new_end\t$new_strand\t$score\n";
          }
        }
      } else {
        # compute best score over all seq_names
        my $seq_regions = {};
        foreach my $seq_name (keys %$filtered) {
          my $mappings = $filtered->{$seq_name};
          my @keys = sort {$mappings->{$a} <=> $mappings->{$b}} keys(%$mappings);
          my $threshold = $mappings->{$keys[0]};
          $seq_regions->{$seq_name} = $threshold;
        }
        my @pot_seq_regions = ();
        my @keys = sort {$seq_regions->{$a} <=> $seq_regions->{$b}} keys(%$seq_regions);
        my $threshold = $seq_regions->{$keys[0]};
        foreach my $seq_region (keys %$seq_regions) {
          if ($seq_regions->{$seq_region} <= $threshold) {
            push @pot_seq_regions, $seq_region;
          }
        }
        foreach my $new_chrom (@pot_seq_regions) {
          my $mappings = $filtered->{$new_chrom};
          my @keys = sort {$mappings->{$a} <=> $mappings->{$b}} keys(%$mappings);
          my $threshold = $mappings->{$keys[0]};
          foreach my $mapping (@keys) {
            my $score = $mappings->{$mapping};
            if ($score <= $threshold) {
              $count_out++;
              my ($new_start, $new_end, $new_strand) = split(':', $mapping);
              push @output_lines, "$id\t$new_chrom\t$new_start\t$new_end\t$new_strand\t$score\n";
            }
          }
        }
      }
    } elsif ($types == 1) {
      
      foreach my $type (keys %{$mapped->{$id}}) {
        my $mappings = {};
        if ($mapped->{$id}->{$type}->{$old_chrom}) {
          foreach my $coords (keys %{$mapped->{$id}->{$type}->{$old_chrom}}) {
            $mappings->{"$old_chrom:$coords"} = $mapped->{$id}->{$type}->{$old_chrom}->{$coords};
          }
        } else {
          foreach my $seq_name (keys %{$mapped->{$id}->{$type}}) {
            foreach my $coords (keys %{$mapped->{$id}->{$type}->{$seq_name}}) {
              $mappings->{"$seq_name:$coords"} = $mapped->{$id}->{$type}->{$seq_name}->{$coords};
            }
          }
        }
        my @keys = sort {$mappings->{$a} <=> $mappings->{$b}} keys(%$mappings);
        my $threshold = $mappings->{$keys[0]};
        foreach my $mapping (@keys) {
          my $score = $mappings->{$mapping};
          if ($score <= $threshold) {
            $count_out++;
            my ($chrom, $new_start, $new_end, $new_strand) = split(':', $mapping);
            push @output_lines, "$id\t$chrom\t$new_start\t$new_end\t$new_strand\t$score\n";
          }
        }
      }
    } else {
      $self->warning("Error for id: $id. More than 2 types");
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
  $self->param('stats_unique_map', $stats_unique_map);
  $self->param('stats_multi_map', $stats_multi_map);
  $self->param('stats_failed', $stats_failed);

}

sub write_statistics {
  my $self = shift;
  my $file_statistics = $self->param('file_statistics');
  my $fh_statistics = FileHandle->new($file_statistics, 'w');

  foreach my $stats (qw/count_input_ids pre_count_unmapped pre_count_mapped stats_failed stats_unique_map stats_multi_map stats_unique_map_after_filter stats_multi_map_after_filter stats_failed_poor_score stats_failed_after_filter stats_failed_non_chrom stats_exceeds_max_map_weight/) {
    if ($self->param($stats)) {
      my $count = $self->param($stats);
      print $fh_statistics "$stats=$count\n";
    }
  }
  $fh_statistics->close();

}

sub join_read_coverage_data {
  my $self = shift;
  my $file_init_feature = $self->param('file_init_feature');
  my $fh_init_feature = FileHandle->new($file_init_feature, 'r'); 

  my $read_coverage_data = {};
  while (<$fh_init_feature>) {
    chomp;
    my $data = $self->read_data($_);
    my $key = $data->{entry};
    $read_coverage_data->{$key} = $data;
  }
  $fh_init_feature->close();

# get new seq_region_ids
  my $seq_region_ids = {};
  my $cdba = $self->param('cdba');
  my $sa = $cdba->get_SliceAdaptor;
  my $slices = $sa->fetch_all('toplevel', undef, 1);
  foreach my $slice (@$slices) {
    $seq_region_ids->{$slice->seq_region_name} = $slice->get_seq_region_id;
  }

# new individual_id
  my $individual_name_oldasm = $self->param('individual_name');
  my $old_individual_id = $self->param('individual_id');

  my $vdba = $self->param('vdba');
  my $ia = $vdba->get_IndividualAdaptor;

  my $individuals_newasm = $ia->fetch_all_by_name($individual_name_oldasm);
#  my $individual_newasm = $individuals_newasm->[0];    
#  my $new_individual_id = $individual_newasm->dbID();
  my $individual_newasm;
  my $new_individual_id;

  if ((scalar @$individuals_newasm) > 1) {
    $individual_newasm = $ia->fetch_by_dbID($old_individual_id); 
    if ($individual_newasm->name eq $individual_name_oldasm) {
      $new_individual_id = $old_individual_id;
    } else {
      die "More than one name for $individual_name_oldasm in new database";
    }
  }
  $individual_newasm = $individuals_newasm->[0];    
  $new_individual_id = $individual_newasm->dbID();

# join feature data with mapping data:
  my $file_load_features = $self->param('file_load_features');
  my $fh_load_features = FileHandle->new($file_load_features, 'w');   
  my $file_filtered_mappings = $self->param('file_filtered_mappings');
  my $fh_mappings = FileHandle->new($file_filtered_mappings, 'r');
  my ($data, $variation_feature_id, $version, $variation_name);
  while (<$fh_mappings>) {
    chomp;
    my ($entry, $seq_name, $start, $end, $strand, $score) = split("\t", $_);
    my $seq_region_id = $seq_region_ids->{$seq_name};
    $data = $read_coverage_data->{$entry};     
    if ($start > $end) {
      $self->warning("Swap start end for $start $end");
      ($start, $end) = ($end, $start);
    }
    $data->{seq_region_id} = $seq_region_id;
    $data->{seq_region_start} = $start;
    $data->{seq_region_end} = $end;
    $data->{individual_id} = $new_individual_id;

    my @output = ();
    foreach my $column_name (sort keys %$data) {
      unless ($column_name =~ /^individual_name$/ || $column_name =~ /^seq_region_name$/ || $column_name =~ /^entry$/) {
        push @output, $data->{$column_name};
      }
    }
    my $line = join("\t", @output);

    print $fh_load_features $line, "\n"; 
  }

  $fh_mappings->close();
  $fh_load_features->close();
}

sub join_svf_data {
  my $self = shift;

  my $file_init_feature = $self->param('file_init_feature');
  my $fh_init_feature = FileHandle->new($file_init_feature, 'r'); 

  my $feature_data = {};
  my $key = '';
  while (<$fh_init_feature>) {
    chomp;
    my $data = $self->read_data($_);
    $key = $data->{structural_variation_feature_id};
    $feature_data->{$key} = $data;
  }
  $fh_init_feature->close();

  # get new seq_region_ids
  my $seq_region_ids = {};
  my $cdba = $self->param('cdba');
  my $sa = $cdba->get_SliceAdaptor;
  my $slices = $sa->fetch_all('toplevel', undef, 1);
  foreach my $slice (@$slices) {
    $seq_region_ids->{$slice->seq_region_name} = $slice->get_seq_region_id;
  }
  # new map_weights
  my $file_filtered_mappings = $self->param('file_filtered_mappings');

  # join feature data with mapping data:
  my $file_load_features = $self->param('file_load_features');
  my $fh_load_features = FileHandle->new($file_load_features, 'w');   
  my $fh_mappings = FileHandle->new($file_filtered_mappings, 'r');
  my ($data, $variation_feature_id, $version, $variation_name);
  while (<$fh_mappings>) {
    chomp;
    my $filtered_data = $self->read_data($_); 
    my $seq_region_name = $filtered_data->{seq_region_name}; 
    my $seq_region_id = $seq_region_ids->{$seq_region_name};
    my $svf_id = $filtered_data->{structural_variation_feature_id}; 
    $data = $feature_data->{$svf_id};
    foreach my $coords (qw/outer_start seq_region_start inner_start inner_end seq_region_end outer_end/) {
      my $location = $data->{coord};
      $data->{seq_region_id} = $seq_region_id || '\N';
    }
    my $line = $self->print_complete_feature_line($data);
    print $fh_load_features $line, "\n"; 
  }

  $fh_mappings->close();
  $fh_load_features->close();
}

sub join_qtl_data {
  my $self = shift;

  my $file_init_feature = $self->param('file_init_feature');
  my $fh_init_feature = FileHandle->new($file_init_feature, 'r'); 

  my $feature_data = {};
  my $key = '';
  while (<$fh_init_feature>) {
    chomp;
    my $data = $self->read_data($_);
    $key = $data->{entry};
    $feature_data->{$key} = $data;
  }
  $fh_init_feature->close();

  # get new seq_region_ids
  my $seq_region_ids = {};
  my $cdba = $self->param('cdba');
  my $sa = $cdba->get_SliceAdaptor;
  my $slices = $sa->fetch_all('toplevel', undef, 1);
  foreach my $slice (@$slices) {
    $seq_region_ids->{$slice->seq_region_name} = $slice->get_seq_region_id;
  }
  # new map_weights
  my $file_filtered_mappings = $self->param('file_filtered_mappings');
  my $fh_mappings = FileHandle->new($file_filtered_mappings, 'r');
  my $map_weights = {};
  while (<$fh_mappings>) {
    chomp;
    my ($query_id, $seq_name, $start, $end, $strand, $score) = split("\t", $_);
    $map_weights->{$query_id}++;
  }
  $fh_mappings->close();

  # join feature data with mapping data:
  my $file_load_features = $self->param('file_load_features');
  my $fh_load_features = FileHandle->new($file_load_features, 'w');   
  $fh_mappings = FileHandle->new($file_filtered_mappings, 'r');
  my ($data, $variation_feature_id, $version, $variation_name);
  while (<$fh_mappings>) {
    chomp;
    my ($query_id, $seq_name, $start, $end, $strand, $score) = split("\t", $_);

    my $seq_region_id = $seq_region_ids->{$seq_name};
    my $map_weight = $map_weights->{$query_id};

    $data = $feature_data->{$query_id};

    $data->{seq_region_id} = $seq_region_id;
    $data->{seq_region_start} = $start;
    $data->{seq_region_end} = $end;
    $data->{seq_region_strand} = $strand;
#    $data->{alignment_quality} = $score;
#    $data->{map_weight} = $map_weight;

    my $line = $self->print_complete_feature_line($data);

    print $fh_load_features $line, "\n"; 
  }

  $fh_mappings->close();
  $fh_load_features->close();
}

sub join_feature_data {
  my $self = shift;

  my $file_init_feature = $self->param('file_init_feature');
  my $fh_init_feature = FileHandle->new($file_init_feature, 'r'); 

  my $feature_data = {};
  my $key = '';
  while (<$fh_init_feature>) {
    chomp;
    my $data = $self->read_data($_);
    if ($self->param('mode') eq 'remap_multi_map') {
      $key = $data->{variation_name};
    } else {
      $key = $data->{variation_feature_id};
    }
    $feature_data->{$key} = $data;
  }
  $fh_init_feature->close();

  # get new seq_region_ids
  my $seq_region_ids = {};
  my $cdba = $self->param('cdba');
  my $sa = $cdba->get_SliceAdaptor;
  my $slices = $sa->fetch_all('toplevel', undef, 1);
  foreach my $slice (@$slices) {
    $seq_region_ids->{$slice->seq_region_name} = $slice->get_seq_region_id;
  }
  # new map_weights
  my $file_filtered_mappings = $self->param('file_filtered_mappings');
  my $fh_mappings = FileHandle->new($file_filtered_mappings, 'r');
  my $map_weights = {};
  while (<$fh_mappings>) {
    chomp;
    my ($query_name, $seq_name, $start, $end, $strand, $score) = split("\t", $_);
    $map_weights->{$query_name}++;
  }
  $fh_mappings->close();

  # join feature data with mapping data:

  my $file_load_features = $self->param('file_load_features');
  my $fh_load_features = FileHandle->new($file_load_features, 'w');   
  $fh_mappings = FileHandle->new($file_filtered_mappings, 'r');
  my ($data, $variation_feature_id, $version, $variation_name);
  while (<$fh_mappings>) {
    chomp;
    my ($query_name, $seq_name, $start, $end, $strand, $score) = split("\t", $_);

    if ($self->param('mode') eq 'remap_multi_map') {
      # query_name: 44919.0-200-1-200-C:C/G:rs2455513
      my @query_name_components = split(':', $query_name);
      $variation_name = pop @query_name_components; 
      $variation_feature_id = shift @query_name_components;
#            my $vf_id_info = shift @query_name_components; 
#            ($variation_feature_id, $version) = split('\.', $vf_id_info);
    } else {
    # query_name: 156358-150-1-150-11:5502587:5502587:1:T/C:rs202026261:dbSNP:SNV
      my @query_name_components = split('-', $query_name, 2);
      $variation_feature_id = $query_name_components[0];
    }
    my $seq_region_id = $seq_region_ids->{$seq_name};
    my $map_weight = $map_weights->{$query_name};

    if ($self->param('mode') eq 'remap_multi_map') {
      $data = $feature_data->{$variation_name};
      $data->{variation_feature_id} = $variation_feature_id;
    } else {
      $data = $feature_data->{$variation_feature_id};
    }

    $data->{seq_region_id} = $seq_region_id;
    $data->{seq_region_start} = $start;
    $data->{seq_region_end} = $end;
    $data->{seq_region_strand} = $strand;
    $data->{alignment_quality} = $score;
    $data->{map_weight} = $map_weight;

    my $line = $self->print_complete_feature_line($data);

    print $fh_load_features $line, "\n"; 
  }

  $fh_mappings->close();
  $fh_load_features->close();
}

sub read_data {
  my $self = shift;
  my $line = shift;
  my @key_values = split("\t", $line);
  my $mapping = {};
  foreach my $key_value (@key_values) {
    my ($table_name, $value) = split('=', $key_value, 2);
    $mapping->{$table_name} = $value;
  }
  return $mapping;
}

sub print_feature_line {
  my ($self, $query_name, $new_seq_info, $score) = @_;
  my ($seq_name, $start, $end, $strand) = split(' ', $new_seq_info);
  my $line =  join("\t", ($query_name, $seq_name, $start, $end, $strand, $score));
  return $line;
}

# returns the row that will be loaded into the new variation_feature_table
# change variation_feature_id to old_variation_feature_id
sub print_complete_feature_line {
  my $self = shift;
  my $data = shift;
  my $variation_feature_id = $data->{variation_feature_id};
  $data->{variation_feature_id_old} = $variation_feature_id;
  my @output = ();
#  $self->warning(join(', ', keys %$data));
  foreach my $column_name (sort keys %$data) {
    unless ($column_name =~ /^variation_feature_id$/ || $column_name =~ /^seq_region_name$/) {
      if ($self->param('mode') eq 'remap_qtls') {
        unless ($column_name =~ /^entry$/ || $column_name =~ /^variation_feature_id_old$/) {
          push @output, $data->{$column_name};
        }
      } else {
        push @output, $data->{$column_name};
      }
    }
  }
  my $line = join("\t", @output);
  return $line;
}

sub get_column_names {
  my ($dbh, $dbname, $table_name) = @_;
  my $query = qq{SHOW columns FROM $dbname.$table_name};
  my $column_names_info = run_query($dbh, $query);
  my @column_names = ();
  foreach my $column_name_info (@$column_names_info) {
    my ($name, $info) = split(',', $column_name_info, 2);
    push @column_names, $name;
  }
  return \@column_names;
}

sub get_table_names {
  my ($dbh, $dbname) = @_;
  my $query = qq{
    SHOW TABLES FROM $dbname;
  };
  my $table_names = run_query($dbh, $query);
  return $table_names;
}

sub run_query {
  my ($dbh, $query) = @_;
  my $sth = $dbh->prepare($query);
  $sth->execute();
  my @results = ();
  while (my $row = $sth->fetchrow_arrayref) {
    my @values = map { defined $_ ? $_ : '\N' } @$row;
    push @results, join(',', @values);
  }
  $sth->finish();
  return \@results;
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

sub get_start {
  my $coord = shift;
  if ($coord) {
    my @keys = keys %$coord;
    return $keys[0];
  }
  return 0;
}

sub overlap {
  my $a = shift;
  my $b = shift;
  return (scalar @$a == scalar @$b);
}

1;
