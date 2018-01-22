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
package Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::PreRunChecks;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::BaseDataDumpsProcess');

use File::Path qw(make_path);

sub run {
  my $self = shift;
  my $file_type  = $self->param('file_type');	
  my $species = $self->param('species');
  my $division = $self->param('division');
  my $pipeline_dir = $self->param_required('pipeline_dir');
  my $script_dir = $self->param_required('script_dir');
  my $tmp_dir = $self->param('tmp_dir');

  # check registry file exists
  # check tmp dir exists

  if ($file_type eq 'gvf') {
    my $gvf_validator = $self->param('gvf_validator');
    my $so_file = $self->param('so_file');
    die "gvf validator not defined" unless (defined $gvf_validator);
    my $gvf_validator_error = `which $gvf_validator 2>&1 1>/dev/null`;
    die "$gvf_validator command not found: $gvf_validator_error" if $gvf_validator_error ne "";
    die "so file not defined" unless (defined $so_file);
    die "wrong location for so file" unless (-f $so_file);
  } elsif ($file_type eq 'vcf') {
    my $vcf_validator = $self->param('vcf_validator');
    my $vcf_validator_error = `which $vcf_validator 2>&1 1>/dev/null`;
    die "$vcf_validator command not found: $vcf_validator_error" if $vcf_validator_error ne "";
    my $vcf_sort = $self->param('vcf_sort');
    my $vcf_sort_error = `which $vcf_sort 2>&1 1>/dev/null`;
    die "$vcf_sort command not found: $vcf_sort_error" if $vcf_sort_error ne "";
  } else {
    die "File type: $file_type is not recognised. It must be gvf or vcf.";
  }

  # Create the tmp_dir
  if (! -d "$tmp_dir") {
    make_path("$tmp_dir") or die "Failed to create dir $tmp_dir $!";
  }

  $self->create_species_dir_tree($species,$division,$pipeline_dir,$script_dir, $file_type);
}

sub write_output {
  my $self = shift;
  $self->dataflow_output_id({}, 2);
  $self->dataflow_output_id({}, 1);
}

sub create_species_dir_tree {
  my ($self,$species,$division, $pipeline_dir, $script_dir, $file_type) = @_;
  my $dump_dir;

  if ($division) {
    my @division = ( ref($division) eq 'ARRAY' ) ? @$division : ($division);
    # If division is defined.
    if ( scalar(@division) ) {
      foreach my $division (@division) {
        $pipeline_dir=$pipeline_dir."/".$division;
        $dump_dir = "$pipeline_dir/$file_type/";
        make_path($dump_dir) unless (-d $dump_dir);
      }
    }
    else{
      $dump_dir = "$pipeline_dir/$file_type/";
      make_path($dump_dir) unless (-d $dump_dir);
    }
  }

  if (-d "$dump_dir/$species") {
    unless (is_empty("$dump_dir/$species")) {
      die("$dump_dir/$species is not empty. Delete files before running the pipeline.");
    }
  } else {
    make_path("$dump_dir/$species") or die "Failed to create dir $dump_dir/$species $!";
  }
}

sub is_empty {
  my $dir = shift;
  opendir(my $dh, $dir) or die "Not a directory $dir";
  my $count =  scalar(grep { $_ ne "." && $_ ne ".." } readdir($dh)) == 0;
  closedir($dh);
  return $count;
}


1;
