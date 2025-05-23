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
package Bio::EnsEMBL::Variation::Pipeline::SpliceAI::FinishFiles;

use strict;
use warnings;
use base ('Bio::EnsEMBL::Variation::Pipeline::SpliceAI::BaseSpliceAI');

use FileHandle;

sub run {
  my $self = shift;
  $self->clean_dir();
}

sub clean_dir {
  my $self = shift;
  my $split_vcf_no_header_dir = $self->param_required('split_vcf_no_header_dir');
  my $split_vcf_input_dir = $self->param_required('split_vcf_input_dir');
  my $split_vcf_output_dir = $self->param_required('split_vcf_output_dir');

  # Remove input files without header
  $self->run_system_command("rm -rf $split_vcf_no_header_dir");

}

1;
