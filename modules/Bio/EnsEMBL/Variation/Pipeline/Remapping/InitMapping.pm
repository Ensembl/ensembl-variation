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
package Bio::EnsEMBL::Variation::Pipeline::Remapping::InitMapping;

use strict;
use warnings;

use FileHandle;
use Bio::DB::Fasta;
use Bio::EnsEMBL::Registry;
use File::Path qw(make_path);

use base ('Bio::EnsEMBL::Variation::Pipeline::Remapping::BaseRemapping');

sub fetch_input {
  my $self = shift;
  my $registry = 'Bio::EnsEMBL::Registry';
  $registry->load_all($self->param('registry_file'));
  my $cdba = $registry->get_DBAdaptor($self->param('species'), 'core');
  my $vdba = $registry->get_DBAdaptor($self->param('species'), 'variation');

  $self->param('registry', $registry);
  $self->param('cdba', $cdba);
  $self->param('vdba', $vdba);
}

sub run {
  my $self = shift;
  my $mode = $self->param('mode');

  if ($mode eq 'remap_multi_map') {
    $self->dump_multi_map_features();
  } 
  elsif ($mode eq 'remap_qtls') {
    $self->dump_features();
    $self->generate_remap_qtls_input();
  }
  elsif ($mode eq 'remap_alt_loci') {
    $self->dump_features_overlapping_alt_loci();
    $self->generate_mapping_input();
  } 
  elsif ($mode eq 'remap_read_coverage') {
    if (!$self->param('use_fasta_files')) {
      $self->dump_read_coverage();
      $self->generate_remap_read_coverage_input();
    }
  } 
  elsif ($mode eq 'remap_svf') {
    $self->dump_features();
    $self->generate_svf_mapping_input();
  } 
  elsif ($mode eq 'remap_post_projection') {
    $self->dump_features();
    $self->generate_mapping_input();
  } 
  else {
    if (!$self->param('use_fasta_files')) {
      $self->dump_features();
      $self->generate_mapping_input();
    } else {
      my $fasta_files_dir = $self->param('fasta_files_dir');
      my $file_count = $self->count_files($fasta_files_dir, 'fa');
      $self->param('file_count', $file_count);
    }
  }
}

sub write_output {
  my $self = shift;
  # initialise mapping jobs
  my $fasta_files_dir = $self->param('fasta_files_dir');
  my $bam_files_dir   = $self->param('bam_files_dir');
  my @jobs = ();

  if ($self->param('mode') eq 'remap_read_coverage') {
    opendir (IND_DIR, $fasta_files_dir) or die $!;
    while (my $individual_dir = readdir(IND_DIR)) {
      next if ($individual_dir =~ /^\./);
      make_path("$bam_files_dir/$individual_dir") or die "Failed to create dir $bam_files_dir/$individual_dir $!";;
      opendir(DIR, "$fasta_files_dir/$individual_dir") or die $!;
      while (my $file = readdir(DIR)) {
        if ($file =~ /^(.+)\.fa$/) {
          my $file_number = $1;
          my $bam_files_dir = "$bam_files_dir/$individual_dir/";
          push @jobs, {
            'file_number'   => $file_number,
             'bam_files_dir' => $bam_files_dir,
             'fasta_file'    => "$fasta_files_dir/$individual_dir/$file_number.fa",
             'sam_file'      => "$bam_files_dir/$file_number.sam",
             'bam_file'      => "$bam_files_dir/$file_number.bam",
             'err_file'      => "$bam_files_dir/$file_number.err",
             'out_file'      => "$bam_files_dir/$file_number.out",
          };
        }
      }
      closedir(DIR);
    }
    closedir(IND_DIR);
  } else {
    my $file_count = $self->param('file_count');
    $self->warning("File count $file_count");
    my $i = 1;
    while ($i <= $file_count) {
      push @jobs, {
        'file_number'   => $i,
          'bam_files_dir' => $bam_files_dir,
          'fasta_file'    => "$fasta_files_dir/$i.fa",
          'sam_file'      => "$bam_files_dir/$i.sam",
          'bam_file'      => "$bam_files_dir/$i.bam",
          'err_file'      => "$bam_files_dir/$i.err",
          'out_file'      => "$bam_files_dir/$i.out",
      };
      $i++;
    }
  }
  $self->dataflow_output_id(\@jobs, 2);
}

sub generate_remap_qtls_input {
  my $self = shift;

  my $old_assembly_fasta_file_dir = $self->param('old_assembly_fasta_file_dir');
  my $fasta_db = Bio::DB::Fasta->new($old_assembly_fasta_file_dir, -reindex => 1);
  $self->param('fasta_db', $fasta_db);

  my $dump_features_dir = $self->param('dump_features_dir');
  my $fasta_files_dir = $self->param('fasta_files_dir');

  my $pipeline_dir = $self->param('pipeline_dir');
  my $fh_report_non_ref_entries = FileHandle->new("$pipeline_dir/report_non_ref_entries.txt", 'w'); 

  my $strand = 1;

  my $file_count = 0;
  opendir(DIR, $dump_features_dir) or die $!;
  while (my $file = readdir(DIR)) {
    if ($file =~ /^(.+)\.txt$/) {
      my $file_number = $1;
      $file_count++;
      my $fh = FileHandle->new("$dump_features_dir/$file", 'r');
      my $fh_fasta_file = FileHandle->new("$fasta_files_dir/$file_number.fa", 'w');
      while (<$fh>) {
        chomp;
        my $data = $self->read_line($_);
        my $seq_region_name = $data->{seq_region_name},
        my $start           = $data->{seq_region_start};
        my $end             = $data->{seq_region_end};
        my $entry           = $data->{entry};
        my $strand = $data->{seq_region_strand};

        unless ($seq_region_name =~ /^\d+$|^X$|^Y$|^MT$/) {
          print $fh_report_non_ref_entries $_, "\n";
        }
        if ($start > $end) {
          $self->warning("End smaller than start for $seq_region_name:$start-$end");
          ($start, $end) = ($end, $start);
        }

        if (($end - $start + 1 ) > 1000) {
          my $upstream_query_sequence   = $self->get_query_sequence($seq_region_name, $start, $start + 500, $strand);
          print $fh_fasta_file ">$entry:upstream\n$upstream_query_sequence\n";
          my $downstream_query_sequence = $self->get_query_sequence($seq_region_name, $end - 500, $end, $strand);
          print $fh_fasta_file ">$entry:downstream\n$downstream_query_sequence\n";
        } elsif (($end - $start) < 100) {
          my $upstream_query_sequence = $self->get_query_sequence($seq_region_name, $start - 100, $start - 1, $strand);
          my $read_sequence = $self->get_query_sequence($seq_region_name, $start, $end, $strand);
          my $downstream_query_sequence = $self->get_query_sequence($seq_region_name, $end + 1, $end + 100, $strand);
          my $read_length = $end - $start + 1; 
          my $id = ">$entry:patched:100:$read_length:100";
          my $sequence = $upstream_query_sequence . $read_sequence . $downstream_query_sequence;    
          print $fh_fasta_file "$id\n$sequence\n";
        } else {
          my $query_sequence = $self->get_query_sequence($seq_region_name, $start, $end, $strand);
          print $fh_fasta_file ">$entry:complete\n$query_sequence\n";
        }
      }
      $fh->close();
      $fh_fasta_file->close();
    }
  }
  $fh_report_non_ref_entries->close();
}

sub generate_remap_read_coverage_input {
  my $self = shift;

  my $old_assembly_fasta_file_dir = $self->param('old_assembly_fasta_file_dir');
  my $fasta_db = Bio::DB::Fasta->new($old_assembly_fasta_file_dir, -reindex => 1);
  $self->param('fasta_db', $fasta_db);

  my $dump_features_dir = $self->param('dump_features_dir');
  my $fasta_files_dir = $self->param('fasta_files_dir');

  my $pipeline_dir = $self->param('pipeline_dir');
  my $fh_report_non_ref_entries = FileHandle->new("$pipeline_dir/report_non_ref_entries.txt", 'w'); 

  my $strand = 1;
  opendir (IND_DIR, $dump_features_dir) or die $!;
  while (my $individual_dir = readdir(IND_DIR)) {
    next if ($individual_dir =~ /^\./);
    make_path("$fasta_files_dir/$individual_dir") or die "Failed to create dir $fasta_files_dir/$individual_dir $!";;
    opendir(DIR, "$dump_features_dir/$individual_dir") or die $!;
    while (my $file = readdir(DIR)) {
      if ($file =~ /^(.+)\.txt$/) {
        my $file_number = $1;
        my $fh = FileHandle->new("$dump_features_dir/$individual_dir/$file", 'r');
        my $fh_fasta_file = FileHandle->new("$fasta_files_dir/$individual_dir/$file_number.fa", 'w');
        while (<$fh>) {
          chomp;
          my $data = $self->read_line($_);
          my $seq_region_name = $data->{seq_region_name},
          my $start           = $data->{seq_region_start};
          my $end             = $data->{seq_region_end};
          my $entry           = $data->{entry};

          unless ($seq_region_name =~ /^\d+$|^X$|^Y$|^MT$/) {
            print $fh_report_non_ref_entries $_, "\n";
          }
          if ($start > $end) {
            $self->warning("End smaller than start for $seq_region_name:$start-$end");
            ($start, $end) = ($end, $start);
          }

          if (($end - $start + 1 ) > 1000) {
            my $upstream_query_sequence   = $self->get_query_sequence($seq_region_name, $start, $start + 500, $strand);
            print $fh_fasta_file ">$entry:upstream\n$upstream_query_sequence\n";
            my $downstream_query_sequence = $self->get_query_sequence($seq_region_name, $end - 500, $end, $strand);
            print $fh_fasta_file ">$entry:downstream\n$downstream_query_sequence\n";
          } elsif (($end - $start) < 100) {
            my $upstream_query_sequence = $self->get_query_sequence($seq_region_name, $start - 100, $start - 1, $strand);
            my $read_sequence = $self->get_query_sequence($seq_region_name, $start, $end, $strand);
            my $downstream_query_sequence = $self->get_query_sequence($seq_region_name, $end + 1, $end + 100, $strand);
            my $read_length = $end - $start + 1; 
            my $id = ">$entry:patched:100:$read_length:100";
            my $sequence = $upstream_query_sequence . $read_sequence . $downstream_query_sequence;    
            print $fh_fasta_file "$id\n$sequence\n";
          } else {
            my $query_sequence = $self->get_query_sequence($seq_region_name, $start, $end, $strand);
            print $fh_fasta_file ">$entry:complete\n$query_sequence\n";
          }
        }
        $fh->close();
        $fh_fasta_file->close();
      }
    }
    closedir(DIR);
  }
  closedir(IND_DIR);
  $fh_report_non_ref_entries->close();
}


sub generate_mapping_input {
  my $self = shift;

  my $old_assembly_fasta_file_dir = $self->param('old_assembly_fasta_file_dir');
  my $fasta_db = Bio::DB::Fasta->new($old_assembly_fasta_file_dir, -reindex => 1);
  $self->param('fasta_db', $fasta_db);

  # store end-coordinates for all seq_regions to check that variation_location + flank_seq_length < slice_end
  my $seq_regions = $self->set_seq_region_boundaries;

  my $dump_features_dir = $self->param('dump_features_dir');
  my $fasta_files_dir   = $self->param('fasta_files_dir');
  my $pipeline_dir      = $self->param('pipeline_dir');
  my $fh_allele_length = FileHandle->new("$pipeline_dir/qc_allele_string_length.txt", 'w');
  my $fh_ref_seq       = FileHandle->new("$pipeline_dir/qc_ref_seq.txt", 'w');

  $self->param('qc_allele_string_length', $fh_allele_length);
  $self->param('qc_ref_seq', $fh_ref_seq);

  my $variants_with_multi_map = {};
  my $mode = $self->param('mode');
  my $dump_multi_map = ($mode eq 'remap_post_projection' || $self->param('dump_multi_map')) ? 1 : 0;
  my $file_count = 0;
  opendir(DIR, $dump_features_dir) or die $!;
  while (my $file = readdir(DIR)) {
    if ($file =~ /^(.+)\.txt$/) {
      my $file_number = $1;
      $file_count++;
      my $fh = FileHandle->new("$dump_features_dir/$file", 'r');
      my $fh_fasta_file = FileHandle->new("$fasta_files_dir/$file_number.fa", 'w');
      while (<$fh>) {
        chomp;
        my $data = $self->read_line($_);
        my $seq_region_name = $data->{seq_region_name},
        my $feature_id      = $data->{variation_feature_id};
        my $start           = $data->{seq_region_start};
        my $end             = $data->{seq_region_end};
        my $strand          = $data->{seq_region_strand};
        my $allele_string   = $data->{allele_string};
        my $map_weight      = $data->{map_weight};
        my $variation_name  = $data->{variation_name};
        if ($map_weight > 1) {
          if ($dump_multi_map) {
            if ($variants_with_multi_map->{$variation_name}) {
              $variants_with_multi_map->{$variation_name}++;
              next;
            } else {
              $variants_with_multi_map->{$variation_name}++;
            }
          } else {
            if ($variants_with_multi_map->{$variation_name}) {
              $variants_with_multi_map->{$variation_name}++;
              next;
            }
          }
        }
        my ($flank_start, $upstream_flank_length, $downstream_flank_length, $flank_end, $variant_length) = @{$self->flank_coordinates($seq_region_name, $start, $end, $strand)};

        # variant surrounded by flank sequences
        my $query_sequence = $self->get_query_sequence($seq_region_name, $flank_start, $flank_end, $strand);

        # replace empty space with underscores <_>
        $allele_string =~ s/\s/_/g;
        $self->qc_alleles($query_sequence, $upstream_flank_length, $variant_length, $allele_string, $variation_name);
        # $allele_string = $self->qc_alleles($config, $query_sequence, $upstream_flank_length, $variant_length, $allele_string, $variation_name);

        my $name = "$seq_region_name:$start:$end:$strand:$variation_name";
        my $id = ">$feature_id-$upstream_flank_length-$variant_length-$downstream_flank_length-$name";
        print $fh_fasta_file "$id\n$query_sequence\n";
      } # end while (read feature file for seq_region)
      $fh->close();
      $fh_fasta_file->close();
    }
  }

  foreach my $file_type (qw/qc_allele_string_length qc_ref_seq/) {
    my $fh = $self->param($file_type);
    $fh->close();
  }

  # store multi map
  my $fh_qc_multi_map = FileHandle->new("$pipeline_dir/multi_map.txt", 'w');

  foreach my $variation_name (keys %$variants_with_multi_map) {
    my $count = $variants_with_multi_map->{$variation_name};
    print $fh_qc_multi_map "$variation_name\t$count\n";
  }
  $fh_qc_multi_map->close();
  $self->param('file_count', $file_count);

}

sub generate_svf_mapping_input {
  my $self = shift;

  my $old_assembly_fasta_file_dir = $self->param('old_assembly_fasta_file_dir');
  my $fasta_db = Bio::DB::Fasta->new($old_assembly_fasta_file_dir, -reindex => 1);
  $self->param('fasta_db', $fasta_db);

  # store end-coordinates for all seq_regions to check that variation_location + flank_seq_length < slice_end
  my $seq_regions = $self->set_seq_region_boundaries;

  my $dump_features_dir = $self->param('dump_features_dir');
  my $fasta_files_dir   = $self->param('fasta_files_dir');
  my $pipeline_dir      = $self->param('pipeline_dir');

  my $file_count = 0;
  opendir(DIR, $dump_features_dir) or die $!;
  while (my $file = readdir(DIR)) {
    if ($file =~ /^(.+)\.txt$/) {
      my $file_number = $1;
      $file_count++;
      my $fh = FileHandle->new("$dump_features_dir/$file", 'r');
      my $fh_lookup = FileHandle->new("$dump_features_dir/lookup_$file_number.txt", 'w');
      my $fh_fasta_file = FileHandle->new("$fasta_files_dir/$file_number.fa", 'w');
      while (<$fh>) {
        chomp;
        my $data = $self->read_line($_);
        my $seq_region_name = $data->{seq_region_name},
        my $feature_id      = $data->{structural_variation_feature_id};
        my $strand          = $data->{seq_region_strand};
        my $variation_name  = $data->{variation_name};
        my @all_coords = ();
        push @all_coords, "seq_region_name=$seq_region_name";
        foreach my $coord_name (qw/outer_start seq_region_start inner_start inner_end seq_region_end outer_end/) {
          my $coord = $data->{$coord_name};
          if ($coord ne '\N') {
            push @all_coords, "$coord_name=$coord";
            my $query_sequence = $self->get_query_sequence($seq_region_name, $coord, $coord + 200, $strand);
            my $id = ">$feature_id-$coord_name";
            print $fh_fasta_file "$id\n$query_sequence\n";
          }
        }
        print $fh_lookup $feature_id, "\t", join(";", @all_coords), "\n";
      } # end while (read feature file for seq_region)
      $fh->close();
      $fh_lookup->close();
      $fh_fasta_file->close();
    }
  }
  $self->param('file_count', $file_count);
}

sub flank_coordinates {
  my ($self, $seq_name, $start, $end, $strand) = @_;

  my $flank_length   = $self->param('flank_seq_length');
  my $flank_start    = 0;
  my $flank_end      = 0;
  my $variant_length = 0;

  my $seq_regions    = $self->param('seq_regions');
  my $slice_end      = $seq_regions->{$seq_name}->{end};

  my $add_to_end = 0;
  my $add_to_start = 0;

  my $upstream_flank_length   = $flank_length;
  my $downstream_flank_length = $flank_length;

  if ($start < $flank_length) {
    $flank_start             = 1;
    $flank_end               = $end + $flank_length + ($flank_length - $start + 1);
    $upstream_flank_length   = $start - 1;
    $downstream_flank_length = $flank_length + ($flank_length - $start + 1);
  } elsif (($end + $flank_length) > $slice_end) {
    my $add_to_start = ($end + $flank_length) - $slice_end;
    $flank_start             = $start - $flank_length - $add_to_start;
    $flank_end               = $slice_end;
    $upstream_flank_length   = $flank_length + $add_to_start;
    $downstream_flank_length = $flank_length - $add_to_start;
  } else {
    $flank_start = $start - $flank_length;
    $flank_end   = $end + $flank_length;
  }

  $variant_length = $end - $start + 1;

  # if insertion
  if ($start > $end) {
    $variant_length = 0;
  }

  return [$flank_start, $upstream_flank_length, $downstream_flank_length, $flank_end, $variant_length];
}

sub get_query_sequence {
  my ($self, $seq_name, $start, $end, $strand) = @_;
  my $fasta_db = $self->param('fasta_db');
  if ($strand == -1) {
    return $fasta_db->seq("$seq_name:$end,$start");
  }
  return $fasta_db->seq("$seq_name:$start,$end");
}

sub qc_alleles {
  my ($self, $query_sequence, $upstream_flank_length, $variant_length, $allele_string, $variation_name) = @_;

  my $upstream   = substr $query_sequence, 0, $upstream_flank_length;
  my $ref        = substr $query_sequence, $upstream_flank_length, $variant_length;
  my $downstream = substr $query_sequence, $upstream_flank_length + $variant_length;

  my @alleles = split('/', $allele_string);
  my $ref_allele = $alleles[0];
  if ($ref_allele eq '-') {
    $ref_allele = '';
  }
  if ($ref ne $ref_allele) {
    my $fh = $self->param('qc_ref_seq');
    print $fh "$variation_name\t$ref_allele\t$ref\n";
  }
}

sub dump_multi_map_features {
  my $self = shift;

  # parse fasta files and store: rsname to fasta_file_number
  # Query_name in fasta_file >14086.1-110-1-497-G:C/G:rs708635
  # parse variation_features on ref_sequence: store variation_feature info in dump_feature/file_number
  # make a note if feature was stored

  my $fasta_files_dir = $self->param('fasta_files_dir');
  my $rs_id_to_file_number = {};
  my $file_numbers = {};
  opendir(DIR, $fasta_files_dir) or die $!;
  while (my $file = readdir(DIR)) {
    if ($file =~ /^(.+)\.fa$/) {
      my $file_number = $1;
      $file_numbers->{$file_number} = 1;
      my $fh = FileHandle->new("$fasta_files_dir/$file", 'r');
      while (<$fh>) {
        chomp;
        my $line = $_;
        if ($line =~ /^>/) {
          my @query_name_components = split(':', $line);
          my $rs_id = pop @query_name_components;
          $rs_id_to_file_number->{$rs_id} = $file_number;
        }
      }
    }
  }
  closedir(DIR);

  my $dump_features_dir = $self->param('dump_features_dir');
  my $fh_hash = {};
  foreach my $file_number (keys %$file_numbers) {
    my $fh = FileHandle->new("$dump_features_dir/$file_number.txt", 'w');
    $fh_hash->{$file_number} = $fh;
  }

  my $file_count = scalar keys %$file_numbers;
  $self->param('file_count', $file_count);

  my $cdba = $self->param('cdba');
  my $vdba = $self->param('vdba');

  my $dbname = $vdba->dbc->dbname();
  my $feature_table = $self->param('feature_table');

  my $dbh = $vdba->dbc->db_handle;
  my $sth = $dbh->prepare(qq{
      SELECT COLUMN_NAME FROM INFORMATION_SCHEMA.COLUMNS
      WHERE TABLE_SCHEMA = '$dbname'
      AND TABLE_NAME = '$feature_table';
      });
  $sth->execute();

  # QC that all necessary columns are there: e.g. seq_region_id, ...
  my @column_names = ();
  while (my @name = $sth->fetchrow_array) {
    push @column_names, $name[0];
  }
  $sth->finish();
  @column_names = sort @column_names;
  my $column_names_concat = join(',', @column_names);
  $self->param('sorted_column_names', $column_names_concat);
  $sth->finish();

  my $sa = $cdba->get_SliceAdaptor;
  # don't include asm exceptions
  my $slices = $sa->fetch_all('toplevel', undef, 0, 1);

  my $seq_region_ids = {};
  foreach my $slice (@$slices) {
    my $seq_region_name = $slice->seq_region_name;
    next if ($seq_region_name =~ /PATCH/);
    my $seq_region_id = $slice->get_seq_region_id;
    $seq_region_ids->{$seq_region_id} = $seq_region_name;
  }

  $sth = $dbh->prepare(qq{
      SELECT variation_name,map_weight,$column_names_concat FROM $feature_table WHERE seq_region_id = ?;
  }, {mysql_use_result => 1});

  my $vf_info_is_stored = {};

  foreach my $seq_region_id (keys %$seq_region_ids) {
    my $seq_region_name = $seq_region_ids->{$seq_region_id};
    $sth->execute($seq_region_id);
    while (my $row = $sth->fetchrow_arrayref) {
      my @values = map { defined $_ ? $_ : '\N' } @$row;
      my $variation_name = shift @values;
      my $map_weight     = shift @values;
      next if ($map_weight == 1);
      next if ($vf_info_is_stored->{$variation_name});
      my $file_number = $rs_id_to_file_number->{$variation_name};
      unless ($file_number) {
        $self->warning("No sequence information for $variation_name from dbSNP");
        next;
      }
      my $fh = $fh_hash->{$file_number};
      my @pairs = ();
      for my $i (0..$#column_names) {
        push @pairs, "$column_names[$i]=$values[$i]";
      }
      push @pairs, "seq_region_name=$seq_region_name";
      print $fh join("\t", @pairs), "\n";
    }
    $sth->finish();
  }

  foreach my $file_number (keys %$fh_hash) {
    my $fh = $fh_hash->{$file_number};
    $fh->close();
  }
}

sub dump_read_coverage {
  my $self = shift;

  my $vdba = $self->param('vdba');
  my $dbname = $vdba->dbc->dbname();
  my $dbh = $vdba->dbc->db_handle;

  # get individiual_ids in old read coverage table
  my $individual_id = 'sample_id';
  my $individual_table = 'sample';
  my $table_names = get_table_names($dbh, $dbname);
  my $rc_individual_ids;
  if (grep /read_coverage/, @$table_names) {
    my $column_names = get_column_names($dbh, $dbname, 'read_coverage');
    if (grep /individual_id/, @$column_names) {
      $individual_id = 'individual_id';
      $individual_table = 'individual';
    }
  # get strain names and number of reads
    my $query = qq{
      SELECT distinct $individual_id FROM read_coverage
    };
    $rc_individual_ids = run_query($dbh, $query);
  }
  my $rc_individual_ids_hash = {};
  foreach my $id (@$rc_individual_ids) {
    $rc_individual_ids_hash->{$id} = 1;
  }

  my $name_to_id = {};

  my @individual_names = split(',', $self->param('individuals'));
  if (scalar @individual_names > 0) {
    foreach my $name (@individual_names) {
      my $ids = $self->get_individual_ids($name);
      my $count = 0;
      foreach my $id (@$ids) {
        if ($rc_individual_ids_hash->{$id}) {
          $name_to_id->{$name} = $id;
          $count++;
        } 
      }
      die "Too many individual ids in table for $name" unless ($count == 1);
    } 
  } else {
    foreach my $id (keys %$rc_individual_ids_hash) {
      my $name = $self->get_individual_name($id);
      $name_to_id->{$name} = $id;
      push @individual_names, $name;
    }
  }

  $individual_id = 'sample_id';
  $individual_table = 'sample';
  $table_names = get_table_names($dbh, $dbname);
  if (grep /read_coverage/, @$table_names) {
    my $column_names = get_column_names($dbh, $dbname, 'read_coverage');
    if (grep /individual_id/, @$column_names) {
      $individual_id = 'individual_id';
      $individual_table = 'individual';
    }
  } else {
    die "No read_coverage table in $dbname";
  }

  my $sth = $dbh->prepare(qq{
      SELECT rc.seq_region_id, rc.seq_region_start, rc.seq_region_end, rc.level, rc.$individual_id, i.name, sr.name
      FROM read_coverage rc, seq_region sr, $individual_table i
      WHERE rc.seq_region_id = sr.seq_region_id
      AND rc.$individual_id = i.$individual_id
      AND i.name = ?; 
      }, {mysql_use_result => 1});
  $sth->execute();

  my @keys = ('seq_region_id', 'seq_region_start', 'seq_region_end', 'level', 'individual_id', 'individual_name', 'seq_region_name');

  my $dump_features_dir = $self->param('dump_features_dir');
  my $entry = 1;
  foreach my $individual_name (@individual_names) {
    my $file_count = 1;
    my $entries_per_file = $self->param('entries_per_file');
    my $count_entries = 0;
    my $individual_id = $name_to_id->{$individual_name};

    unless (-d "$dump_features_dir/$individual_id") {
      make_path("$dump_features_dir/$individual_id") or die "Failed to create dir $dump_features_dir/$individual_id $!";;
    }

    my $fh = FileHandle->new("$dump_features_dir/$individual_id/$file_count.txt", 'w');

    $sth->execute($individual_name);
    while (my $row = $sth->fetchrow_arrayref) {
      my @values = map { defined $_ ? $_ : '\N' } @$row;
      my @pairs = ();
      for my $i (0..$#keys) {
        push @pairs, "$keys[$i]=$values[$i]";
      }
      push @pairs, "entry=$entry";
      $entry++;
      if ($count_entries >= $entries_per_file) {
        $fh->close();
        $file_count++;
        $fh = FileHandle->new("$dump_features_dir/$individual_id/$file_count.txt", 'w');
        $count_entries = 0;
      }
      $count_entries++;
      print $fh join("\t", @pairs), "\n";
    }
    $fh->close();
    $sth->finish();
  }
}

sub dump_features_overlapping_alt_loci {
  my $self = shift; 

  my $cdba = $self->param('cdba');

  my $alt_loci_to_coords = {};
  my $ref_to_alt_loci = {};
  my $ref_to_unique_region_coords = {};

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
          $alt_loci_to_coords->{$seq_region_name}->{start} = $start;
          $alt_loci_to_coords->{$seq_region_name}->{end} = $end;
        }
      } elsif ($assembly_exception_type eq 'REF') {
        my $assembly_exception_features = $aefa->fetch_all_by_Slice($slice);
        my $ref_start = $slice->start;
        my $ref_end = $slice->end; 
        $ref_to_unique_region_coords->{$seq_region_name}->{start} = $ref_end;
        $ref_to_unique_region_coords->{$seq_region_name}->{end} = $ref_start;
        foreach my $feature (@$assembly_exception_features) {
          my $alt_slice = $feature->alternate_slice();
          my $alt_slice_name = $alt_slice->seq_region_name;
          unless ($alt_slice_name =~ /^X$|^Y$|PATCH/) {
            $ref_to_alt_loci->{$seq_region_name}->{$alt_slice_name} = 1; 
          }
        }
      }
  }
# compute boundaries
  foreach my $ref (keys %$ref_to_alt_loci) {
    foreach my $alt_loci (keys %{$ref_to_alt_loci->{$ref}}) {
      my $start = $alt_loci_to_coords->{$alt_loci}->{start};
      my $end = $alt_loci_to_coords->{$alt_loci}->{end};
      if ($ref_to_unique_region_coords->{$ref}->{start} > $start) {
        $ref_to_unique_region_coords->{$ref}->{start} = $start;
      }
      if ($ref_to_unique_region_coords->{$ref}->{end} < $end) {
        $ref_to_unique_region_coords->{$ref}->{end} = $end;
      }
    }
  } 

  my $vdba = $self->param('vdba');

  my $dbname = $vdba->dbc->dbname();
  my $feature_table = $self->param('feature_table');

  my $dbh = $vdba->dbc->db_handle;
  my $sth = $dbh->prepare(qq{
      SELECT COLUMN_NAME FROM INFORMATION_SCHEMA.COLUMNS
      WHERE TABLE_SCHEMA = '$dbname'
      AND TABLE_NAME = '$feature_table';
      });
  $sth->execute();

# QC that all necessary columns are there: e.g. seq_region_id, ...
  my @column_names = ();
  while (my @name = $sth->fetchrow_array) {
    push @column_names, $name[0];
  }
  $sth->finish();
  @column_names = sort @column_names;
  my $column_names_concat = join(',', @column_names);
  $self->param('sorted_column_names', $column_names_concat);
  $sth->finish();

  my $dump_features_dir = $self->param('dump_features_dir');
# don't include asm exceptions
  $slices = $sa->fetch_all('toplevel', undef, 0, 1);

  my $seq_region_ids = {};
  foreach my $slice (@$slices) {
    my $seq_region_name = $slice->seq_region_name;
    if ($slice->assembly_exception_type eq 'REF') {
      my $seq_region_id = $slice->get_seq_region_id;
      $seq_region_ids->{$seq_region_id} = $seq_region_name;
    }
  }

  $sth = $dbh->prepare(qq{
      SELECT seq_region_start,seq_region_end,source_id,$column_names_concat FROM $feature_table WHERE seq_region_id = ?;
      }, {mysql_use_result => 1});

  my $file_count = 1;
  my $entries_per_file = $self->param('entries_per_file');
  my $count_entries = 0;
  my $fh = FileHandle->new("$dump_features_dir/$file_count.txt", 'w');

  foreach my $seq_region_id (keys %$seq_region_ids) {
    my $seq_region_name = $seq_region_ids->{$seq_region_id};
    my @alt_loci_names = keys %{$ref_to_alt_loci->{$seq_region_name}};
    $sth->execute($seq_region_id);
    while (my $row = $sth->fetchrow_arrayref) {
      my @values = map { defined $_ ? $_ : '\N' } @$row;
      my $seq_region_start = shift @values;
      my $seq_region_end = shift @values;
      my $source_id = shift @values;

      next unless ($source_id == 1);

      my $overlaps_alt_loci = 0;
      foreach my $alt_loci_name (@alt_loci_names) {
        my $alt_loci_start = $alt_loci_to_coords->{$alt_loci_name}->{start};
        my $alt_loci_end   = $alt_loci_to_coords->{$alt_loci_name}->{end};
        if ($seq_region_start >= $alt_loci_start && $seq_region_end <= $alt_loci_end) {
          $overlaps_alt_loci = 1;
          last;
        } 
      }
#        my $alt_loci_start = $ref_to_unique_region_coords->{$seq_region_name}->{start};
#        my $alt_loci_end = $ref_to_unique_region_coords->{$seq_region_name}->{end};
#        next unless ($seq_region_start >= $alt_loci_start && $seq_region_end <= $alt_loci_end);
      next unless ($overlaps_alt_loci);
      my @pairs = ();
      for my $i (0..$#column_names) {
        push @pairs, "$column_names[$i]=$values[$i]";
      }
      push @pairs, "seq_region_name=$seq_region_name";
      if ($count_entries >= $entries_per_file) {
        $fh->close();
        $file_count++;
        $fh = FileHandle->new("$dump_features_dir/$file_count.txt", 'w');
        $count_entries = 0;
      }
      $count_entries++;
      print $fh join("\t", @pairs), "\n";
    }
    $sth->finish();
  }
  $fh->close();
  $self->param('file_count', $file_count);

}

sub dump_features {
  my $self = shift;
  my $cdba = $self->param('cdba');
  my $vdba = $self->param('vdba');

  my $feature_table = $self->param('feature_table');

  my $dbname = $vdba->dbc->dbname();
  my $dbh = $vdba->dbc->db_handle;
  my $sth = $dbh->prepare(qq{
      SELECT COLUMN_NAME FROM INFORMATION_SCHEMA.COLUMNS
      WHERE TABLE_SCHEMA = '$dbname'
      AND TABLE_NAME = '$feature_table';
      });
  $sth->execute();

# QC that all necessary columns are there: e.g. seq_region_id, ...
  my @column_names = ();
  while (my @name = $sth->fetchrow_array) {
    push @column_names, $name[0];
  }
  $sth->finish();
  @column_names = sort @column_names;
  my $column_names_concat = join(',', @column_names);
  $self->param('sorted_column_names', $column_names_concat);
  $sth->finish();

  my $dump_features_dir = $self->param('dump_features_dir');
  my $sa = $cdba->get_SliceAdaptor;
# don't include asm exceptions
  my $slices = $sa->fetch_all('toplevel', undef, 0, 1);

  my $seq_region_ids = {};
  foreach my $slice (@$slices) {
    my $seq_region_name = $slice->seq_region_name;
    next if ($seq_region_name =~ /PATCH/);
    my $seq_region_id = $slice->get_seq_region_id;
    $seq_region_ids->{$seq_region_id} = $seq_region_name;
  }

  my $file_count = 1;
  my $entries_per_file = $self->param('entries_per_file');
  my $count_entries = 0;
  my $fh = FileHandle->new("$dump_features_dir/$file_count.txt", 'w');

  my $mode = $self->param('mode');
  my @tables = ('feature_table');
  if ($self->param('feature_table_failed_projection')) {
    push @tables, 'feature_table_failed_projection'
  }


  my $entry = 1;

  foreach my $table (@tables) {
    my $feature_table = $self->param($table); 
    if (($mode eq 'remap_post_projection') && ($table eq 'feature_table')) {
      $sth = $dbh->prepare(qq{
        SELECT $column_names_concat FROM $feature_table WHERE seq_region_id = ? AND map_weight > 1;
      }, {mysql_use_result => 1});
    } else {
      $sth = $dbh->prepare(qq{
        SELECT $column_names_concat FROM $feature_table WHERE seq_region_id = ?;
      }, {mysql_use_result => 1});
    }
    foreach my $seq_region_id (keys %$seq_region_ids) {
      my $seq_region_name = $seq_region_ids->{$seq_region_id};
      $sth->execute($seq_region_id) or die $sth->errstr;
      while (my $row = $sth->fetchrow_arrayref) {
        my @values = map { defined $_ ? $_ : '\N' } @$row;
        my @pairs = ();
        for my $i (0..$#column_names) {
          push @pairs, "$column_names[$i]=$values[$i]";
        }
        push @pairs, "seq_region_name=$seq_region_name";
        if ($feature_table eq 'phenotype_feature') {
          push @pairs, "entry=$entry";
          $entry++;
        }

        if ($count_entries >= $entries_per_file) {
          $fh->close();
          $file_count++;
          $fh = FileHandle->new("$dump_features_dir/$file_count.txt", 'w');
          $count_entries = 0;
        }
        $count_entries++;
        print $fh join("\t", @pairs), "\n";
      }
      $sth->finish();
    }
  }
  $fh->close();
  $self->param('file_count', $file_count);
}

sub get_individual_ids {
  my $self = shift;
  my $individual_name = shift;
  $self->warning($individual_name);
  my $vdba = $self->param('vdba');

  my $dbname = $vdba->dbc->dbname();
  my $dbh = $vdba->dbc->db_handle();

  my $column = 'sample_id';
  my $table = 'sample';

  my $column_names = get_column_names($dbh, $dbname, 'individual');
  if (grep /individual_id/, @$column_names) {
    $column = 'individual_id';
    $table = 'individual';
  }
  my $sth = $dbh->prepare(qq{SELECT $column FROM $table where name='$individual_name';});

  my @ids = ();
  $sth->execute();
  while (my $row = $sth->fetchrow_arrayref) {
    push @ids, $row->[0];
  }
  $sth->finish();
  return \@ids;
}

sub get_individual_name {
  my $self = shift;
  my $individual_id = shift;
  $self->warning($individual_id);
  my $vdba = $self->param('vdba');

  my $dbname = $vdba->dbc->dbname();
  my $dbh = $vdba->dbc->db_handle();

  my $column = 'sample_id';
  my $table = 'sample';

  my $column_names = get_column_names($dbh, $dbname, 'individual');
  if (grep /individual_id/, @$column_names) {
    $column = 'individual_id';
    $table = 'individual';
  }
  my $sth = $dbh->prepare(qq{SELECT name FROM $table where $column=$individual_id;});

  my @names = ();
  $sth->execute();
  while (my $row = $sth->fetchrow_arrayref) {
    push @names, $row->[0];
  }
  $sth->finish();
  return $names[0];
}

sub set_seq_region_boundaries {
  my $self = shift;
  my $cdba = $self->param('cdba');
  my $sa = $cdba->get_SliceAdaptor;
  # don't include asm exceptions, fetch the full length of the Y chromosome
  my $slices = $sa->fetch_all('toplevel', undef, 0, 1);
  my $seq_regions = {};
  foreach my $slice (@$slices) {
    my $end             = $slice->end;
    my $seq_region_name = $slice->seq_region_name;
    if ($end > $seq_regions->{$seq_region_name}->{end}) {
      $seq_regions->{$seq_region_name}->{end} = $end;
    }
  }
  $self->param('seq_regions', $seq_regions);
  return $seq_regions;
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


1;
