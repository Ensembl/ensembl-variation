=head1 LICENSE
Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2025] EMBL-European Bioinformatics Institute
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

package Bio::EnsEMBL::Variation::Pipeline::Remapping::InitReadCoverageMapping;

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
  $self->dump_read_coverage();
  $self->generate_remap_read_coverage_input();
}

sub write_output {
  my $self = shift;
  # initialise mapping jobs
  my $fasta_files_dir = $self->param('fasta_files_dir');
  my $bam_files_dir   = $self->param('bam_files_dir');
  my @jobs = ();

  opendir (SAMPLE_DIR, $fasta_files_dir) or die $!;
  while (my $sample_dir = readdir(SAMPLE_DIR)) {
    next if ($sample_dir =~ /^\./);
    make_path("$bam_files_dir/$sample_dir") or die "Failed to create dir $bam_files_dir/$sample_dir $!";;
    opendir(DIR, "$fasta_files_dir/$sample_dir") or die $!;
    while (my $file = readdir(DIR)) {
      if ($file =~ /^(.+)\.fa$/) {
        my $file_number = $1;
        my $bam_files_dir = "$bam_files_dir/$sample_dir/";
        push @jobs, {
          'file_number'   => $file_number,
          'bam_files_dir' => $bam_files_dir,
          'fasta_file'    => "$fasta_files_dir/$sample_dir/$file_number.fa",
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
  $self->dataflow_output_id(\@jobs, 2);
}


sub generate_remap_read_coverage_input {
  my $self = shift;

  my $old_assembly_fasta_file_dir = $self->param('old_assembly_dir');
  my $fasta_db = Bio::DB::Fasta->new($old_assembly_fasta_file_dir, -reindex => 1);
  $self->param('fasta_db', $fasta_db);

  my $dump_features_dir = $self->param('dump_features_dir');
  my $fasta_files_dir = $self->param('fasta_files_dir');

  my $pipeline_dir = $self->param('pipeline_dir');
  my $fh_report_non_ref_entries = FileHandle->new("$pipeline_dir/report_non_ref_entries.txt", 'w'); 

  my $strand = 1;
  opendir (SAMPLE_DIR, $dump_features_dir) or die $!;
  while (my $sample_dir = readdir(SAMPLE_DIR)) {
    next if ($sample_dir =~ /^\./);
    make_path("$fasta_files_dir/$sample_dir") or die "Failed to create dir $fasta_files_dir/$sample_dir $!";;
    opendir(DIR, "$dump_features_dir/$sample_dir") or die $!;
    while (my $file = readdir(DIR)) {
      if ($file =~ /^(.+)\.txt$/) {
        my $file_number = $1;
        my $fh = FileHandle->new("$dump_features_dir/$sample_dir/$file", 'r');
        my $fh_fasta_file = FileHandle->new("$fasta_files_dir/$sample_dir/$file_number.fa", 'w');
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

sub dump_read_coverage {
  my $self = shift;

  my $vdba = $self->param('vdba_oldasm');
  my $dbname = $vdba->dbc->dbname();
  my $dbh = $vdba->dbc->db_handle;

  my $query = qq{
    SELECT distinct sample_id FROM read_coverage
  };
  my $sample_ids = $self->run_query($dbh, $query);

  my $sth = $dbh->prepare(qq{
      SELECT rc.seq_region_id, rc.seq_region_start, rc.seq_region_end, rc.level, rc.sample_id, s.name, sr.name
      FROM read_coverage rc, seq_region sr, sample s
      WHERE rc.seq_region_id = sr.seq_region_id
      AND rc.sample_id = s.sample_id
      AND s.sample_id = ?;
      }, {mysql_use_result => 1});
  $sth->execute();

  my @keys = ('seq_region_id', 'seq_region_start', 'seq_region_end', 'level', 'sample_id', 'sample_name', 'seq_region_name');

  my $dump_features_dir = $self->param('dump_features_dir');
  my $entry = 1;
  foreach my $sample_id (@$sample_ids) {
    my $file_count = 1;
    my $entries_per_file = $self->param('entries_per_file');
    my $count_entries = 0;

    unless (-d "$dump_features_dir/$sample_id") {
      make_path("$dump_features_dir/$sample_id") or die "Failed to create dir $dump_features_dir/$sample_id $!";;
    }

    my $fh = FileHandle->new("$dump_features_dir/$sample_id/$file_count.txt", 'w');

    $sth->execute($sample_id);
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
        $fh = FileHandle->new("$dump_features_dir/$sample_id/$file_count.txt", 'w');
        $count_entries = 0;
      }
      $count_entries++;
      print $fh join("\t", @pairs), "\n";
    }
    $fh->close();
    $sth->finish();
  }
}

1;
