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

package Bio::EnsEMBL::Variation::Pipeline::Remapping::RunMapping;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Hive::Process');

use FileHandle;

sub fetch_input {
    my $self = shift;
    1;
}

sub run {
    my $self = shift;

    my $tool_dir      = $self->param('tool_dir'); 	
    my $bwa_dir       = $self->param('bwa_dir');
    my $samtools_dir  = $self->param('samtools_dir');
    my $fasta_file    = $self->param('fasta_file');
    my $sam_file      = $self->param('sam_file');
    my $file_number   = $self->param('file_number');
    my $bam_files_dir = $self->param('bam_files_dir');

    my $sam2bam_err = "$bam_files_dir/sam2bam.$file_number.err";
    my $sam2bam_out = "$bam_files_dir/sam2bam.$file_number.out";

    my $bam_file = $self->param('bam_file');
    my $err_file = $self->param('err_file');
    my $out_file = $self->param('out_file');

    my $new_assembly_fasta_file_dir = $self->param('new_assembly_fasta_file_dir');
    my $new_assembly_fasta_file_name = $self->param('new_assembly_fasta_file_name');
    my $look_up_file = "$new_assembly_fasta_file_dir/$new_assembly_fasta_file_name";

#my $cmd = "$tool_dir/$bwa_dir/bwa mem -a $look_up_file $fasta_file";
#    my $cmd = "$tool_dir/$bwa_dir/bwa mem -a $look_up_file $fasta_file";
    my $cmd = "/homes/anja/bin/bwa/bwa mem -a $look_up_file $fasta_file";


    $self->run_cmd("$cmd 1>$sam_file 2>$err_file");

    $cmd = "gzip $sam_file";
    $self->run_cmd($cmd);

#    $cmd = "$tool_dir/$samtools_dir/samtools view -uS $sam_file.gz | $tool_dir/$samtools_dir/samtools sort - $bam_files_dir/$file_number"; 
    $cmd = "$samtools_dir view -uS $sam_file.gz | $samtools_dir sort -o $bam_files_dir/$file_number.bam -"; 

#    $cmd = "$samtools_dir view -uS $sam_file.gz | $samtools_dir sort sort -T /tmp/aln.sorted -o aln.sorted.bam aln.bam       - $bam_files_dir/$file_number"; 

#  $cmd = "$samtools_dir sort -O bam -l 0 -o $bam_files_dir/$file_number $sam_file.gz";



    $self->run_cmd("$cmd 1>$sam2bam_out 2>$sam2bam_err");

    return 1;
}

sub run_cmd {
    my $self = shift;
    my $cmd = shift;
    if (my $return_value = system($cmd)) {
        $return_value >>= 8;
        die "system($cmd) failed: $return_value";
    }
}

1;
