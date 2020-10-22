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
package Bio::EnsEMBL::Variation::Pipeline::SpliceAI::MergeFiles;

use strict;
use warnings;
use base ('Bio::EnsEMBL::Variation::Pipeline::SpliceAI::BaseSpliceAI');

use FileHandle;

sub run {
  my $self = shift;
  $self->merge_vcf_files();
}

sub merge_vcf_files {
  my $self = shift;
  my $input_dir = $self->param_required('input_dir');
  my $chr_dir = $self->param_required('chr_dir');
  my $output_dir = $self->param_required('output_dir');
  my $output_file_name = $self->param_required('output_file_name');

  my $input_dir_chr = $input_dir . '/' . $chr_dir . '/vcf_files';
  $self->param('input_dir_chr', $input_dir_chr);

  opendir(my $read_dir, $input_dir_chr) or die $!;

  while(my $tmp_vcf = readdir($read_dir)) {
    next if ($tmp_vcf =~ m/^\./ || $tmp_vcf =~ m/\.gz/);

    # Some variants don't have a score - the main reason is because the transcript is not in the gene annotation file
    # Before merging the files, the variants without scores need to be deleted
    # $self->run_system_command("sed -i \"/\t\.\t\.\t\./d\" $input_dir_chr/$tmp_vcf");

    $self->run_system_command("bgzip $input_dir_chr/$tmp_vcf");
    $self->run_system_command("tabix -p vcf $input_dir_chr/$tmp_vcf.gz");

  }
  close($read_dir);

  # $self->run_system_command("ls $input_dir_chr/*.vcf.gz | split -l 10 - $input_dir_chr/split2merge");
  # $self->run_system_command("for f in $input_dir_chr/split2merge*; do bcftools merge -l \$f -m none -Oz -o tmp_\$f.vcf.gz; done");
  # $self->run_system_command("bcftools merge -m none -Oz -o $output_file_name$chr_dir.vcf.gz $input_dir_chr/tmp_*.vcf.gz");
}

1;
