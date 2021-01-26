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
package Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::InitJoinDump;

use strict;
use warnings;

use FileHandle;

use base ('Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::BaseDataDumpsProcess');

sub run {
  my $self = shift;
  my $species      = $self->param('species');
  my $pipeline_dir = $self->data_dir($species);
  my $file_type = $self->param('file_type');
  my $working_dir = "$pipeline_dir/$file_type/$species/";
  my $mode = $self->param('mode');
  if ($mode eq 'join_split_slice') {
    my $input = $self->get_split_slice_files($working_dir, $file_type);
    $self->param('input', $input);
  } else { 
    my @input = ();
    foreach my $file_type (qw/gvf vcf/) {

      my $dir = "$pipeline_dir/$file_type/$species/";
      my $files = {};

      opendir(my $dh, $dir) or die $!;
      my @dir_content = readdir($dh);
      closedir($dh);
      foreach my $file (@dir_content) {
        next if ($file =~ m/^\./);
        # homo_sapiens_structural_variations-chr10.gvf
        # homo_sapiens_structural_variations-chr10.vcf.gz 
        if ($file =~ m/\.($file_type\.gz|$file_type)$/) {
          my $file_name = $file;
          $file_name =~ s/\.($file_type\.gz|$file_type)$//;
          # $file_name e.g. homo_sapiens_structural_variations-chr10
          my ($dump_type, $range) = split('-', $file_name);
          push @{$files->{$dump_type}}, $range if ($range);
        }
      }
      foreach my $dump_type (keys %$files) {
        my $file_name = $dump_type;
        next if ($dump_type eq 'homo_sapiens_generic' || $dump_type eq 'homo_sapiens_incl_consequences');
        if (scalar @{$files->{$dump_type}} > 0) {
          push @input, {
            dir => $dir,
            file_name => $file_name,
            dump_type => $dump_type,
            input_ids => $files->{$dump_type},
            file_type => $file_type,
            mode => 'final_join',
          };
        }
      }
    }
    $self->param('input', \@input);
  }
}

# Get file components for slice_split files which look like
# homo_sapiens_generic-131555_146405131_151285302.gvf
# file_name-seq_region_id-start-end.gvf
sub get_split_slice_files {
  my ($self, $working_dir, $file_type) = @_;
  my $dump_types = {};
  my @input = ();
  opendir(my $dh, $working_dir) or die $!;
  my @dir_content = readdir($dh);
  closedir($dh);
  foreach my $file (@dir_content) {
    if ($file =~ m/([a-z|_]+)-([0-9]+)_([0-9]+)_([0-9]+)\.$file_type/) {
      my $dump_type = $1;
      my $seq_region_id = $2;
      $dump_types->{$dump_type}->{$seq_region_id} = 1;
    }
  }
  foreach my $dump_type (keys %$dump_types) {
    foreach my $seq_region_id (keys %{$dump_types->{$dump_type}}) {
      push @input, {
        mode => 'join_split_slice',
        file_type => $file_type,
        working_dir => $working_dir,
        dump_type => $dump_type,
        seq_region_id => $seq_region_id,
      };
    }
  }
  return \@input;
}

sub write_output {
  my $self = shift;
  if (scalar @{$self->param('input')} > 0) {
    $self->dataflow_output_id($self->param('input'), 2);
  } else {
    $self->dataflow_output_id([{mode => 'no_join'}], 2);
  }
  return;
}

1;
