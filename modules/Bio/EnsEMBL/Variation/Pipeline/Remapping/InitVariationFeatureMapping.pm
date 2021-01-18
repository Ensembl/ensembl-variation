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

package Bio::EnsEMBL::Variation::Pipeline::Remapping::InitVariationFeatureMapping;

use strict;
use warnings;

use FileHandle;
use Bio::DB::Fasta;
use Bio::EnsEMBL::Registry;
use File::Path qw(make_path);

use base ('Bio::EnsEMBL::Variation::Pipeline::Remapping::InitMapping');

sub fetch_input {
  my $self = shift;
  $self->SUPER::fetch_input;
}

sub run {
  my $self = shift;
  $self->SUPER::dump_features();
  $self->generate_mapping_input();
}

sub write_output {
  my $self = shift;
  $self->SUPER::write_output;
}

sub generate_mapping_input {
  my $self = shift;

  my $old_assembly_fasta_file_dir = $self->param('old_assembly_dir');
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
  my $dump_multi_map = 1;
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
        if (!$query_sequence) {
          $self->warning("no query sequence for $seq_region_name, $flank_start, $flank_end");
        }

        # replace empty space with underscores <_>
        $allele_string =~ s/\s/_/g;
        $self->qc_alleles($query_sequence, $upstream_flank_length, $variant_length, $allele_string, $variation_name);
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

1;
