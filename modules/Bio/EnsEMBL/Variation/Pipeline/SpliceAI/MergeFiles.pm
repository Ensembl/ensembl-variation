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

  # Final file + final file after it is sorted
  my $final_file = $output_dir . '/' . $output_file_name . $chr_dir . '.vcf';
  my $final_file_sorted = $output_dir . '/sorted_' . $output_file_name . $chr_dir . '.vcf';

  # Write final file header
  # Header contains the INFO for SpliceAI scores
  open(my $write, '>', $final_file) or die $!;
  my $header_line = "##fileformat=VCFv4.2\n##contig=<ID=1,length=248956422>\n##contig=<ID=2,length=242193529>\n".
    "##contig=<ID=3,length=198295559>\n##contig=<ID=4,length=190214555>\n##contig=<ID=5,length=181538259>\n".
    "##contig=<ID=6,length=170805979>\n##contig=<ID=7,length=159345973>\n##contig=<ID=8,length=145138636>\n".
    "##contig=<ID=9,length=138394717>\n##contig=<ID=10,length=133797422>\n##contig=<ID=11,length=135086622>\n".
    "##contig=<ID=12,length=133275309>\n##contig=<ID=13,length=114364328>\n##contig=<ID=14,length=107043718>\n".
    "##contig=<ID=15,length=101991189>\n##contig=<ID=16,length=90338345>\n##contig=<ID=17,length=83257441>\n".
    "##contig=<ID=18,length=80373285>\n##contig=<ID=19,length=58617616>\n##contig=<ID=20,length=64444167>\n##contig=<ID=21,length=46709983>\n".
    "##contig=<ID=22,length=50818468>\n##contig=<ID=X,length=156040895>\n##contig=<ID=Y,length=57227415>\n##contig=<ID=MT,length=16569>\n".
    "##INFO=<ID=SpliceAI,Number=.,Type=String,Description=\"SpliceAIv1.3 variant annotation. These include delta scores (DS) and delta positions (DP) for acceptor gain (AG), acceptor loss (AL), donor gain (DG), and donor loss (DL). Format: ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL\">\n".
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
  print $write $header_line;
  close($write);

  opendir(my $read_dir, $input_dir_chr) or die $!;

  while(my $tmp_vcf = readdir($read_dir)) {
    next if ($tmp_vcf =~ m/^\./);

    # Append the file content (excluding the headers) into final file
    $self->run_system_command("grep \"^[^#]\" $input_dir_chr/$tmp_vcf >> $final_file ");

  }
  close($read_dir);

  # Sort final file
  $self->run_system_command("bcftools sort -o $final_file_sorted $final_file");
}

1;
