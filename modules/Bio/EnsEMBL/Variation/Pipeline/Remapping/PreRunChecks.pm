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
package Bio::EnsEMBL::Variation::Pipeline::Remapping::PreRunChecks;

use strict;
use warnings;

use Bio::EnsEMBL::Registry;
use File::Path qw(make_path remove_tree);
use FileHandle;
use Bio::EnsEMBL::Registry;
use IPC::Cmd qw(can_run);

use base ('Bio::EnsEMBL::Variation::Pipeline::Remapping::BaseRemapping');

sub fetch_input {
  my $self = shift;
}

sub run {
  my $self = shift;

  my $pipeline_dir = $self->param('pipeline_dir');
  die "$pipeline_dir doesn't exist" unless (-d $pipeline_dir);		

  foreach my $folder (qw/old_assembly_dir new_assembly_dir/) {
    my $dir = $self->param($folder);
    die "$dir for $folder doesn't exist" unless (-d "$dir");		
  }

  my $dir = $self->param('new_assembly_dir');
  foreach my $file_type (('.fa.amb', '.fa.ann', '.fa.bwt', '.fa.pac', '.fa.sa')) {
    unless ($self->count_files($dir, $file_type)) {
      die("New assembly file is not indexed. $file_type is missing.");
    }
  }

  foreach my $tool (qw/bwa samtools/) {
    can_run($self->param($tool)) or die "$tool could not be found at location: ", $self->param($tool);
  }

  foreach my $file (qw/ensembl_regsitry_oldasm ensembl_regsitry_newasm/) {
    "File $file is missing " if (! -f $self->param($file));
  }

  my @folders = qw/bam_files_dir filtered_mappings_dir load_features_dir mapping_results_dir statistics_dir dump_mapped_features_dir/;
  push @folders, (qw/qc_failure_reasons_dir qc_mapped_features_dir qc_update_features_dir/) if ($self->param('mode') eq 'remap_variation_feature');
  foreach my $folder (@folders) {
    my $dir = $self->param($folder);
    if (-d $dir) {
      remove_tree($dir);
    }
    make_path($dir);
  } 
  
  if (!$self->param('use_fasta_files')) {
    foreach my $folder (qw/dump_features_dir fasta_files_dir/) {
      my $dir = $self->param($folder);
      if (-d $dir) {
        remove_tree($dir);
      }
      make_path($dir);
    } 
  } else {
    my $dir = $self->param('fasta_files_dir');
    my $count = $self->count_files($dir, '.fa');
    if ($count == 0) {
      die ("There are no fasta_files. Set parameter 'generate_fasta_files' to 1 in the conf file.");
    }				
  # remove index files
    if ($self->param('mode') eq 'remap_read_coverage') {
      opendir (IND_DIR, $dir) or die $!;
      while (my $individual_dir = readdir(IND_DIR)) {
        next if ($individual_dir =~ /^\./);
        $self->run_cmd("rm $dir/$individual_dir/*.fai");
      }
      closedir (IND_DIR);
    } else {
      $self->run_cmd("rm -f $dir/*.fai");
    }
  }
}

sub write_output {
  my $self = shift;
}

1;
