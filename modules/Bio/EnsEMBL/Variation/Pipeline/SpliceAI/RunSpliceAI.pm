=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2020] EMBL-European Bioinformatics Institute

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
package Bio::EnsEMBL::Variation::Pipeline::SpliceAI::RunSpliceAI;

use strict;
use warnings;
use base ('Bio::EnsEMBL::Variation::Pipeline::SpliceAI::BaseSpliceAI');

use FileHandle;
use Bio::EnsEMBL::IO::Parser::BedTabix;
use Bio::EnsEMBL::IO::Parser::VCF4Tabix;

sub run {
  my $self = shift;
  $self->set_chr_from_filename();
  $self->run_spliceai();
}

sub run_spliceai {
  my $self = shift;
  my $main_dir = $self->param_required('main_dir');
  my $vcf_input_dir = $self->param('split_vcf_input_dir'); # $self->o('main_dir') . '/tmp_split_vcf_input'
  my $vcf_file = $self->param('input_file');
  my $output_dir = $self->param_required('output_dir');
  my $fasta_file = $self->param_required('fasta_file');
  my $gene_annotation = $self->param_required('gene_annotation');

  my $chr = $self->param('chr');

  my $vcf_input_dir_chr = $vcf_input_dir.'/chr'.$chr;

  if (! -d $vcf_input_dir_chr) {
    die("Directory ($vcf_input_dir_chr) doesn't exist");
  }

  my $output_dir_chr = $output_dir."/chr".$chr;
  my $out_files_dir = $output_dir_chr."/out_files";
  my $output_vcf_files_dir = $output_dir_chr."/vcf_files";

  if (! -d $output_vcf_files_dir) {
    die("Directory ($output_vcf_files_dir) doesn't exist");
  }

  my $err = $out_files_dir."/".$vcf_file.".err";
  my $out = $out_files_dir."/".$vcf_file.".out";

  my $cmd = "spliceai -I $vcf_input_dir_chr/$vcf_file -O $output_vcf_files_dir/$vcf_file -R $fasta_file -A $gene_annotation";
  $self->run_system_command($cmd);
  
  $self->warning('Error: ' . $stderr . ' Code: ' . $exit_code);
}

1;
