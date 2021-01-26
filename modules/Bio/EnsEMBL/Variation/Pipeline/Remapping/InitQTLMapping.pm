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

package Bio::EnsEMBL::Variation::Pipeline::Remapping::InitQTLMapping;

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
  $self->dump_qtls();
  $self->generate_remap_qtls_input();
}

sub write_output {
  my $self = shift;
  $self->SUPER::write_output;
}

sub generate_remap_qtls_input {
  my $self = shift;

  my $old_assembly_fasta_file_dir = $self->param('old_assembly_dir');
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
        my $entry           = $data->{phenotype_feature_id};
        my $strand          = $data->{seq_region_strand};

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

sub dump_qtls {
  my $self = shift;
  $self->SUPER::dump_features("AND type = 'QTL'");
}

1;
