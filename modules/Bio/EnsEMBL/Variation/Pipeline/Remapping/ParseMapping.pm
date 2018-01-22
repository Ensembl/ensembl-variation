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
package Bio::EnsEMBL::Variation::Pipeline::Remapping::ParseMapping;

use strict;
use warnings;

#use Bio::DB::Sam;
use Bio::DB::HTS;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Registry;
use FileHandle;

use base ('Bio::EnsEMBL::Hive::Process');


sub fetch_input {
  my $self = shift;
  my $file_number         = $self->param('file_number');
  my $bam_file            = $self->param('bam_file');
  my $fasta_file          = $self->param('fasta_file');
  my $registry_file       = $self->param('registry_file');
  my $species             = $self->param('species');
  my $mapping_results_dir = $self->param('mapping_results_dir');
  my $compare_locations  = $self->param('compare_locations');

  my $sam = Bio::DB::HTS->new( -bam => $bam_file, -fasta => $fasta_file,);	
  $self->param('sam', $sam);	

  if ($self->param('mode') eq 'remap_read_coverage') {
    my $individual_id = $self->param('individual_id');
    my $fh_mappings        = FileHandle->new("$mapping_results_dir/$individual_id/mappings_$file_number.txt", 'w');
    my $fh_failed_mappings = FileHandle->new("$mapping_results_dir/$individual_id/failed_mapping_$file_number.txt", 'w');

    $self->param('fh_mappings', $fh_mappings);
    $self->param('fh_failed_mappings', $fh_failed_mappings);
  } else {
    my $fh_mappings        = FileHandle->new("$mapping_results_dir/mappings_$file_number.txt", 'w');
    my $fh_failed_mappings = FileHandle->new("$mapping_results_dir/failed_mapping_$file_number.txt", 'w');

    $self->param('fh_mappings', $fh_mappings);
    $self->param('fh_failed_mappings', $fh_failed_mappings);
  }
  if ($compare_locations) {
    my $registry = 'Bio::EnsEMBL::Registry';
    $registry->load_all($registry_file);
    my $vfa = $registry->get_adaptor($species, 'variation', 'variationfeature'); 
    $self->param('vfa', $vfa);
  }
}

sub run {
  my $self = shift;
  my $mode = $self->param('mode');
  if ($mode eq 'remap_read_coverage' || $mode eq 'remap_svf' || $mode eq 'remap_qtls') {
    $self->parse_read_location();
  } else {
    $self->parse_variation_location();
  }
}

sub parse_read_location {
  my $self = shift;
  my $sam                = $self->param('sam');
  my $fh_mappings        = $self->param('fh_mappings');
  my $fh_failed_mappings = $self->param('fh_failed_mappings');
  my $fasta_file         = $self->param('fasta_file');

  my @alignments = $sam->features();

  my $map_weights = $self->get_map_weights(\@alignments);
  $self->test_all_variants_are_mapped($fasta_file, $map_weights);
  my $query_sequences = $self->get_query_sequences(\@alignments, $map_weights);

  foreach my $alignment (@alignments) {
    my $query_name = $alignment->query->name;
    my $seq_region_name = $alignment->seq_id;
    if ($seq_region_name =~ m/:/) {
      my @seq_id_parts = split(':', $seq_region_name);
      $seq_region_name = $seq_id_parts[2];
    }

    my $t_start = $alignment->start;
    my $t_end   = $alignment->end;
    my $q_start = $alignment->query->start;
    my $q_end   = $alignment->query->end;

    my $t_strand = $alignment->strand;
    my $diff = $t_end - $t_start + 1;
    my $map_weight = $map_weights->{$query_name};		

    my @components = split(':', $query_name);
    if ((scalar @components) == 5) { # query sequence is patch
      my $id = $components[0];
      my $type = $components[1];
      my $upstream_seq_length = $components[2];
      my $read_length = $components[2];
      $t_start = $t_start + 100 - 1;
      $t_end = $t_start + $read_length - 1;
      $query_name = "$id:$type";
    }
    unless ($seq_region_name && $t_start && $t_end){
      print $fh_failed_mappings join("\t", 'NO_MAPPING'), "\n";
      next;
    } 

    my @cigar_string = @{$alignment->cigar_array};

    my $clipped_nucleotides = 0;
    foreach my $sub_string (@cigar_string) {
      my $operation = $sub_string->[0];
      my $count     = $sub_string->[1];
      if ($operation eq 'H' || $operation eq 'S') {
        $clipped_nucleotides += $count;
      }	
    }

    my $cigar_string = $alignment->cigar_str;
    my $edit_distance         = $alignment->aux_get("NM");
    my $alignment_score_count = $alignment->aux_get("AS"); 

    my $query_seq = $alignment->query->dna; 	
    if ($query_sequences->{$query_name}) {
      $query_seq = $query_sequences->{$query_name};
    }

    my $length_query_seq      = length($query_seq);

    my $relative_alignment_score = 0;
    if ($length_query_seq == 0) {
      print $fh_failed_mappings join("\t", 'NO_MAPPING', $query_name, $seq_region_name, $t_start, $t_end, $q_start, $q_end, $t_strand, $map_weight, 0, $diff, $cigar_string), "\n";
      next;
    } else {
      $relative_alignment_score = ($length_query_seq - ($clipped_nucleotides + $edit_distance)) / $length_query_seq;	
      $relative_alignment_score = sprintf("%.5f", $relative_alignment_score);
    }
    if ($relative_alignment_score > 0.8) {
      print $fh_mappings join("\t", $query_name, $seq_region_name, $t_start, $t_end, $t_strand, $map_weight, $relative_alignment_score, $cigar_string), "\n";
    } else {
      print $fh_failed_mappings join("\t", 'LOW_SCORE', $query_name, $seq_region_name, $t_start, $t_end, $q_start, $q_end, $t_strand, $map_weight, $relative_alignment_score, $diff, $cigar_string), "\n";
    }

  }
  $fh_mappings->close();
  $fh_failed_mappings->close();

}

sub parse_variation_location {
  my $self = shift;
  $self->warning("Parse variation location");
  my $sam                = $self->param('sam');
  my $fh_mappings        = $self->param('fh_mappings');
  my $fh_failed_mappings = $self->param('fh_failed_mappings');
  my $vfa                = $self->param('vfa');
  my $fasta_file         = $self->param('fasta_file');
  my $compare_locations  = $self->param('compare_locations');

  my ($variation_name, $vf_id, $length_before, $length_var, $length_after);	

  my @alignments = $sam->features();

  my $map_weights = $self->get_map_weights(\@alignments);
  $self->test_all_variants_are_mapped($fasta_file, $map_weights);
  my $query_sequences = $self->get_query_sequences(\@alignments, $map_weights);

  foreach my $alignment (@alignments) {
    my $query_name = $alignment->query->name;
    my ($vf_id, $length_before_var, $length_var, $length_after, $variation_name) = split('-', $query_name, 5);
    my $q_strand = $alignment->strand;
    if ($q_strand != 1) {
      $length_before_var = $length_after;
    }	

    my $map_weight = $map_weights->{$query_name};		
    my @cigar_string = @{$alignment->cigar_array};

    my $clipped_nucleotides = 0;
    foreach my $sub_string (@cigar_string) {
      my $operation = $sub_string->[0];
      my $count     = $sub_string->[1];
      if ($operation eq 'H' || $operation eq 'S') {
        $clipped_nucleotides += $count;
      }	
    }

    my $query_seq = $alignment->query->dna; 	
    if ($query_sequences->{$query_name}) {
      $query_seq = $query_sequences->{$query_name};
    }
    my $results = $self->map_variant($alignment, $length_before_var, $length_var);
    my ($snp_t_start, $snp_t_end, $indel_colocation, $flag_suspicious) = @$results;
    unless ($snp_t_start && $snp_t_end && $q_strand){
      $flag_suspicious = 1;
    } 

    my $vf;
    # old seq info
    my $old_seq_info = '';
    if ($compare_locations) {	
      eval {$vf = $vfa->fetch_by_dbID($vf_id);};
      if ($@) {
        die "Error fetching VF $vf_id for $query_name): $@";
      }
      my $seq_region_name = $vf->seq_region_name;
      my $vf_start        = $vf->seq_region_start;
      my $vf_end          = $vf->seq_region_end;
      my $strand          = $vf->strand;
      $old_seq_info = join(" ", ($seq_region_name, $vf_start, $vf_end, $strand));
    }
    # new seq info	
    my $new_seq_id = $alignment->seq_id;
    if ($new_seq_id =~ m/:/) {
      my @seq_id_parts = split(':', $new_seq_id);
      $new_seq_id = $seq_id_parts[2];
    }
    my $new_seq_info = join(" ", ($new_seq_id, $snp_t_start, $snp_t_end, $q_strand));

    my $cigar                 = $alignment->cigar_str;
    my $edit_distance         = $alignment->aux_get("NM");
    my $alignment_score_count = $alignment->aux_get("AS"); 

    my $length_query_seq      = length($query_seq);
    my $count_ns              = $self->count_number_of_ns_in_clipped_seq($alignment, $query_seq);


    my $relative_alignment_score = 0;
    if ($length_query_seq == 0) {
      $flag_suspicious = 1;
    } else {
      $relative_alignment_score = ($length_query_seq - ($clipped_nucleotides + $edit_distance)) / $length_query_seq;	
      $relative_alignment_score = sprintf("%.5f", $relative_alignment_score);
    }

    my $mapping_info = join("\t", (
          $query_name, $map_weight, $cigar, $relative_alignment_score, "clipped nucleotides $clipped_nucleotides"));
    my $mapping_result = join("\t", ($old_seq_info, $new_seq_info, $mapping_info));
    if ($relative_alignment_score > 0.5) { # should be 0.5
      unless ($flag_suspicious || $indel_colocation) {
        print $fh_mappings $mapping_result, "\n";
      } else {
        my $pre = $indel_colocation ? "INDEL\t" : "SUSPICIOUS\t";
        print $fh_failed_mappings $pre, $mapping_result, "\n";
      }
    } else {
      print $fh_failed_mappings "LOW_ALGN_SCORE\t", $mapping_result, "\n";
    }
  }

  $fh_mappings->close();
  $fh_failed_mappings->close();
}

sub map_variant {
  my $self = shift;
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
    # warn
    #        next;
    # return
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
    } elsif ($operation eq 'I')	{
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

sub test_all_variants_are_mapped {
  my $self = shift;
  my ($fasta_file, $map_weights) = @_;
  # test same number ids as in input fasta file
  my $input_ids = {};
  my $fh_fasta_file = FileHandle->new($fasta_file, 'r');
  while (<$fh_fasta_file>) {
    chomp;
    if (/^>/) {
      my $query_name = $_;
      $query_name =~ s/>//;
      $input_ids->{$query_name} = 1;
    }
  }
  $fh_fasta_file->close();

  my $count_mapped_ids = scalar keys %$map_weights;
  my $count_input_ids = scalar keys %$input_ids;

  unless ($count_mapped_ids == $count_input_ids) {
    my $unmapped_input_ids = {};
    foreach my $id (keys %$map_weights) {
      unless ($input_ids->{$id}) {
        $self->warning("Not mapped $id");
      }
    }
    die("Number of ids differ in fasta ($count_input_ids) and bam file ($count_mapped_ids): $fasta_file");
  }
}

sub get_map_weights {
  my $self = shift;
  my $alignments = shift;
  my $map_weights;
  foreach my $alignment (@$alignments) {
    $map_weights->{$alignment->query->name}++;
  }
  return $map_weights;
}

sub count_number_of_ns_in_clipped_seq {
  my $self = shift;
  my $alignment = shift;
  my $query_seq = shift;
  my $query_seq_length = length($query_seq);
  my @cigar_string = @{$alignment->cigar_array};
  my $pos_in_query_seq = 1; 
  my $count_ns = 0;
  while (my $sub_string = shift @cigar_string) {
    my $operation = $sub_string->[0];
    my $count     = $sub_string->[1];
    if ($operation eq 'H' || $operation eq 'S') {
      my $substr = substr $query_seq, $pos_in_query_seq - 1, $count;
      if (length($substr) == $count)	{
        my @count_ns = split(//, $substr);
        foreach my $letter (@count_ns) {
          if ($letter eq 'N') {
            $count_ns++;
          }			
        }
      } else {
        my $debug_cigar_str = $alignment->cigar_str;
        my $debug_query_name = $alignment->query->name;
        die("Length of substring differs from count in cigar string for $query_seq $debug_cigar_str $debug_query_name");
      }
      $pos_in_query_seq += $count;
    } elsif ($operation eq 'M') {
      $pos_in_query_seq += $count;
    } elsif ($operation eq 'I') {
      $pos_in_query_seq += $count;                    
    } elsif ($operation eq 'D') {
    # no change in position
    }
  }
  return $count_ns; 
}

sub get_query_sequences {
  my $self = shift;
  my $alignments = shift;
  my $map_weights = shift;
  my $query_sequences;	
  # to save space BWA abbreviates sequences for query sequences that map multiple times:
  foreach my $alignment (@$alignments) {
    my $query_name = $alignment->query->name;
    if ($map_weights->{$query_name} > 1 ) {
      my $query_sequence = $alignment->query->dna; 
      unless ($query_sequences->{$query_name}) {
        $query_sequences->{$query_name} = $query_sequence;
      } else {
        if (length($query_sequences->{$query_name}) < length($query_sequence)) {
          $query_sequences->{$query_name} = $query_sequence;
        }
      }
    }
  }   
  return $query_sequences;
}

sub write_output {
  my $self = shift;
  $self->dataflow_output_id({
      'fh_mappings'        => $self->param('fh_mappings'),
      'fh_failed_mappings' => $self->param('fh_failed_mappings'),
  }, 1);	
}

1;
