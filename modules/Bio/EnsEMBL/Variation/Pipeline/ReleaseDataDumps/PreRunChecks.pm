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
package Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::PreRunChecks;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::BaseDataDumpsProcess');

use File::Path qw(make_path);

sub run {
  my $self = shift;
  my $division = $self->param('division');
  my $pipeline_dir = $self->param_required('pipeline_dir');
  my $script_dir = $self->param_required('script_dir');
  my $registry_file = $self->param_required('ensembl_registry');
  my $tmp_dir = $self->param('tmp_dir');

  die "pipeline_dir ($pipeline_dir) doesn't exist" unless (-d $pipeline_dir);
  die "script_dir ($script_dir) doesn't exist" unless (-d $script_dir);
  die "Registry file ($registry_file) doesn't exist" unless (-f $registry_file);

  if (! -d "$tmp_dir") {
    make_path("$tmp_dir") or die "Failed to create dir $tmp_dir $!";
  }

  my $gvf_validator = $self->param('gvf_validator');
  die "gvf validator not defined" unless (defined $gvf_validator);
  my $gvf_validator_error = `which $gvf_validator 2>&1 1>/dev/null`;
  die "$gvf_validator command not found: $gvf_validator_error" if $gvf_validator_error ne "";

  my $so_file = $self->param('so_file');
  die "so file not defined" unless (defined $so_file);
  die "wrong location for so file" unless (-f $so_file);

  my $vcf_validator = $self->param('vcf_validator');
  my $vcf_validator_error = `which $vcf_validator 2>&1 1>/dev/null`;
  die "$vcf_validator command not found: $vcf_validator_error" if $vcf_validator_error ne "";

  my $vcf_sort = $self->param('vcf_sort');
  my $vcf_sort_error = `which $vcf_sort 2>&1 1>/dev/null`;
  die "$vcf_sort command not found: $vcf_sort_error" if $vcf_sort_error ne "";

  $self->create_species_dir_tree($division, $pipeline_dir);
}

sub write_output {
  my $self = shift;
}

sub create_species_dir_tree {
  my ($self, $division, $pipeline_dir) = @_;
  if ($division) {
    my @division = ( ref($division) eq 'ARRAY' ) ? @$division : ($division);
    foreach my $division (@division) {
      $pipeline_dir=$pipeline_dir."/".$division."/variation";
      make_path($pipeline_dir) unless (-d $pipeline_dir);
    }
  }
}

1;
