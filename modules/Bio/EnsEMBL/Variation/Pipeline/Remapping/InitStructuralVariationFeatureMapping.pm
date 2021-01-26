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

package Bio::EnsEMBL::Variation::Pipeline::Remapping::InitStructuralVariationFeatureMapping;

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
  $self->generate_svf_mapping_input();
}

sub write_output {
  my $self = shift;
  $self->SUPER::write_output;
}

sub generate_svf_mapping_input {
  my $self = shift;

  my $old_assembly_fasta_file_dir = $self->param('old_assembly_dir');
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
          if ($coord ne '\N') { # some coord types (outer_*, inner_*) are not always defined for structural variations
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

1;
