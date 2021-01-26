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

=head1 NAME

Bio::EnsEMBL::Variation::Pipeline::DBSNPImport::LoadDBSNPFile

=head1 DESCRIPTION

Loads a dbSNP JSON file

=cut

package Bio::EnsEMBL::Variation::Pipeline::DBSNPImport::LoadDBSNPFile;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);

use FileHandle;
use Bio::EnsEMBL::Registry;
use File::Path qw(make_path);
use Bio::EnsEMBL::Variation::Utils::Date;
use POSIX;

sub fetch_input {
  my $self = shift;
}

sub run {
  my $self = shift;
  
  # The filename is flowing in as refsnp-chr21-aa.gz
  # The input_directory is <data-dir>/<sub_dir>
  
  my $filename = $self->param_required('filename');
  my $sub_dir = $self->param_required('sub_dir');
  my $data_dir = $self->param_required('data_dir');
  my $rpt_dir = $self->param_required('rpt_dir');
  my $registry_file = $self->param_required('ensembl_registry');
  my $ancestral_fasta_file = $self->param_required('ancestral_fasta_file');
  my $fasta_file = $self->param_required('fasta_file');
  my $script_dir = $self->param_required('script_dir');
  my $assembly = $self->param_required('assembly');

  $self->warning("filename ($filename)");
  $self->warning("subdir ($sub_dir)");
  $self->warning("data_dir ($data_dir)");
  $self->warning("rpt_dir ($rpt_dir)");
  $self->warning("registry_file ($registry_file)");
  $self->warning("ancestral_fasta_file ($ancestral_fasta_file)");
  $self->warning("fasta_file ($fasta_file)");
  $self->warning("script_dir ($script_dir)");
  $self->warning("assembly ($assembly)");
  
  my $load_script = join("/", $script_dir, 'load_dbsnp.pl');
  my $data_dir_run = join("/", $data_dir, $sub_dir);
  my $rpt_dir_run = join("/", $rpt_dir, $sub_dir);

  $self->warning("data_dir_run = $data_dir_run");
  $self->warning("rpt_dir_run = $rpt_dir_run");
  $self->warning("load_script = $load_script");
 
  my $input_file = join("/", $data_dir_run, $filename);
  my $cmd;

  $cmd = join(" ", 'perl',  $load_script ,
            '-registry' ,   $registry_file ,
            '-input_file', $input_file,
            '-rpt_dir' ,   $rpt_dir_run,
            '-ancestral_fasta_file' , $ancestral_fasta_file,
            '-fasta_file', $fasta_file,
            '-assembly', $assembly);
  
  $self->warning($cmd);
  # Do the system call
  my ($return_value, $stderr, $flat_cmd) = $self->run_system_command($cmd);
  if ($return_value) {
      die("there was an error running as ($flat_cmd: $stderr)");
  }  
}

1;
