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
=cut

package Bio::EnsEMBL::Variation::Pipeline::Remapping::InitMapping;

use strict;
use warnings;

use FileHandle;
use Bio::EnsEMBL::Registry;
use File::Path qw(make_path);

use base ('Bio::EnsEMBL::Variation::Pipeline::Remapping::BaseRemapping');

sub fetch_input {
  my $self = shift;
  my $registry = 'Bio::EnsEMBL::Registry';
  $registry->load_all($self->param('registry_file_oldasm'));
  my $cdba = $registry->get_DBAdaptor($self->param('species'), 'core');
  my $vdba = $registry->get_DBAdaptor($self->param('species'), 'variation');

  $self->param('registry_oldasm', $registry);
  $self->param('cdba_oldasm', $cdba);
  $self->param('vdba_oldasm', $vdba);

  if ($self->param('seq_region_name_mappings_file')) {
    my $fh = FileHandle->new($self->param('seq_region_name_mappings_file'), 'r');
    my $seq_region_name_mappings = {};
    while (<$fh>) {
      chomp;
      my ($fasta_id, $ensembl_id) = split/\s/;
      $seq_region_name_mappings->{$ensembl_id} = $fasta_id;
    }
    $self->param('seq_region_name_mappings', $seq_region_name_mappings);
  }
 $self->test_bioperl_version;
}

sub run {
  my $self = shift;
}

sub dump_features {
  my $self = shift;
  my $extra_sql = shift; # AND type = 'QTL'
  my $cdba_oldasm = $self->get_oldasm_core_database_connection;
  my $vdba_oldasm = $self->get_oldasm_variation_database_connection;
  my $feature_table = $self->param('feature_table');

  my @column_names = @{$self->get_sorted_column_names($vdba_oldasm, $feature_table)};
  my $column_names_string = join(',', @column_names);
  $self->param('sorted_column_names', $column_names_string);

  my $dump_features_dir = $self->param('dump_features_dir');
  my $file_count = 1;
  my $entries_per_file = $self->param('entries_per_file');
  my $count_entries = 0;
  my $fh = FileHandle->new("$dump_features_dir/$file_count.txt", 'w');

  my $dbh = $vdba_oldasm->dbc->db_handle();
  my $sth = $dbh->prepare(qq{
    SELECT $column_names_string FROM $feature_table WHERE seq_region_id = ? $extra_sql;
  }, {mysql_use_result => 1});

  my $seq_region_ids = $self->get_seq_region_ids($cdba_oldasm);

  foreach my $seq_region_id (keys %$seq_region_ids) {
    my $seq_region_name = $seq_region_ids->{$seq_region_id};
    if ($self->param('debug')) {
      next if ("$seq_region_name" ne $self->param('debug_sequence_name'));
    }
    $sth->execute($seq_region_id) or die $sth->errstr;
    while (my $row = $sth->fetchrow_arrayref) {
      my @values = map { defined $_ ? $_ : '\N' } @$row;
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

sub write_output {
  my $self = shift;
  # initialise mapping jobs
  my $fasta_files_dir = $self->param('fasta_files_dir');
  my $bam_files_dir   = $self->param('bam_files_dir');
  my @jobs = ();

  my $file_count = $self->param('file_count');
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
  $self->dataflow_output_id(\@jobs, 2);
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
  my $seq_region_name_mappings = $self->param('seq_region_name_mappings');  
  if ($seq_region_name_mappings) {
    $seq_name = $seq_region_name_mappings->{$seq_name};
  }
  return undef if (!$seq_name);
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

sub set_seq_region_boundaries {
  my $self = shift;
  my $cdba = $self->param('cdba_oldasm');
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

1;
