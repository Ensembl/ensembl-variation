=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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
package Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::FileUtils;

use strict;

use base ('Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::BaseDataDumpsProcess');

sub fetch_input {
  my $self = shift;
}

sub run {
  my $self = shift;
  my $mode = $self->param('mode');
  $self->post_gvf_dump_cleanup if ($mode eq 'post_gvf_dump_cleanup');
  $self->post_gvf2vcf_cleanup if ($mode eq 'post_gvf2vcf_cleanup');
}

sub post_gvf_dump_cleanup {
  my $self = shift;
  my $data_dump_dir = $self->param('pipeline_dir');
  my $tmp_dir = $self->param('tmp_dir');
  my $species = $self->param('species');

  system("gzip $data_dump_dir/gvf/$species/*.gvf");
  system("cat $data_dump_dir/gvf/$species/Validate_* > $tmp_dir/GVF_Validate_$species");
  system("rm $data_dump_dir/gvf/$species/Validate_*");
  system("cat $data_dump_dir/gvf/$species/*.{err,out} > $tmp_dir/GVF_$species");
  system("rm $data_dump_dir/gvf/$species/*.{err,out,txt}");
}

sub post_gvf2vcf_cleanup {
  my $self = shift;
  my $data_dump_dir = $self->param('pipeline_dir');
  my $tmp_dir = $self->param('tmp_dir');
  my $species = $self->param('species');
  system("cat $data_dump_dir/vcf/$species/Validate_* > $tmp_dir/VCF_Validate_$species");
  system("rm $data_dump_dir/vcf/$species/Validate_*");
  system("cat $data_dump_dir/vcf/$species/*.{err,out} > $tmp_dir/VCF_$species");
  system("rm $data_dump_dir/vcf/$species/*.{err,out}");
  system("rm $data_dump_dir/vcf/$species/*_validate.vcf.gz");
  system("rm $data_dump_dir/vcf/$species/*.vcf");
}

sub write_output {
  my $self = shift;
}

1;
