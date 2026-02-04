=head1 LICENSE
Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2026] EMBL-European Bioinformatics Institute
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
package Bio::EnsEMBL::Variation::Pipeline::SpliceAI::SpliceAI_conf;

use strict;
use warnings;

use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
 # All Hive databases configuration files should inherit from HiveGeneric, directly or indirectly
use base ('Bio::EnsEMBL::Hive::PipeConfig::EnsemblGeneric_conf');

sub default_options {
    my ($self) = @_;

    # The hash returned from this function is used to configure the
    # pipeline, you can supply any of these options on the command
    # line to override these default values.
    
    # You shouldn't need to edit anything in this file other than
    # these values, if you find you do need to then we should probably
    # make it an option here, contact the variation team to discuss
    # this - patches are welcome!

    return {
        %{ $self->SUPER::default_options()
        },    # inherit other stuff from the base class

        hive_auto_rebalance_semaphores => 1,
        hive_default_max_retry_count => 0,
        hive_force_init      => 1,
        hive_use_param_stack => 1,

        ensembl_release            => $self->o('ensembl_release'),
        assembly                   => $self->o('assembly'),
        pipeline_name              => 'spliceai_scores',
        main_dir                   => $self->o('main_dir'), # main directory where all files and directories are going to be stored
        input_directory            => $self->o('main_dir') . '/input_vcf_files', # input files
        split_vcf_no_header_dir    => $self->o('main_dir') . '/split_vcf_no_header', # contains the input files after being splitted (files without headers)
        split_vcf_input_dir        => $self->o('main_dir') . '/split_vcf_input', # contains the splitted input vcf files with headers, these are the files used to run SpliceAI
        split_vcf_output_dir       => $self->o('main_dir') . '/split_vcf_output', # temporary output files, still splitted
        output_dir                 => $self->o('main_dir') . '/output', # final output files already merged by chromosome
        fasta_file                 => $self->o('fasta_file'),
        gene_annotation            => $self->o('gene_annotation'),
        step_size                  => 2000000, # number of variants used to split the main vcf files
        check_transcripts          => 0, # if set to 1 checks which are the new MANE Select transcripts for the last months and only calculates SpliceAI scores for these variants overlapping these transcripts
        transcripts_from_file      => undef,
        time_interval              => 4, # checks which transcripts were updated/created in the last 4 months; only used if check_transcripts = 1 and we want to check the new transcripts in the core db
        masked_scores              => 1, # calculate masked scores
        registry                   => undef, # database where new MANE transcripts are going to be checked; only used if check_transcripts = 1
        output_file_name           => 'spliceai_final_scores_' . $self->o('ensembl_release') . '_',

        pipeline_wide_analysis_capacity => 80,

        pipeline_db => {
            -host   => $self->o('hive_db_host'),
            -port   => $self->o('hive_db_port'),
            -user   => $self->o('hive_db_user'),
            -pass   => $self->o('hive_db_password'),
            -dbname => $ENV{'USER'} . '_ehive_' . $self->o('pipeline_name') . '_' . $self->o('ensembl_release') . '_' . $self->o('assembly'),
            -driver => 'mysql',
        },
    };
}

sub resource_classes {
    my ($self) = @_;
    return {
        %{$self->SUPER::resource_classes},
        'gpu'      => {
                        'SLURM' => '--time=24:00:00 --gres=gpu:1 --mem=32G'
                      },
        '4Gb_job'  => {
                        'SLURM' => "--partition=standard --time=4:00:00 --mem=8G"
                      },
         'default' => {
                        'SLURM' => "--partition=standard --time=1:00:00 --mem=4G"
                      }
    };
}

sub pipeline_analyses {
  my ($self) = @_;
  my @analyses;
  push @analyses, (
      {   -logic_name => 'init_files',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
          -input_ids  => [{}],
          -parameters => {
            'input_directory' => $self->o('input_directory'),
            'inputcmd'        => 'find #input_directory# -type f -name "all_snps_ensembl_*.vcf.gz" -printf "%f\n"',
          },
          -flow_into  => {
            '2->A' => {'split_files' => {'vcf_file' => '#_0#'}},
            'A->1' => ['get_chr_dir']
          },
      },
      { -logic_name => 'split_files',
        -module => 'Bio::EnsEMBL::Variation::Pipeline::SpliceAI::SplitFiles',
        -input_ids  => [],
        -rc_name => '4Gb_job',
        -parameters => {
          'main_dir'                   => $self->o('main_dir'),
          'input_directory'            => $self->o('input_directory'),
          'split_vcf_no_header_dir'    => $self->o('split_vcf_no_header_dir'),
          'split_vcf_input_dir'        => $self->o('split_vcf_input_dir'),
          'split_vcf_output_dir'       => $self->o('split_vcf_output_dir'),
          'step_size'                  => $self->o('step_size'),
          'check_transcripts'          => $self->o('check_transcripts'),
          'registry'                   => $self->o('registry'),
          'transcripts_from_file'      => $self->o('transcripts_from_file'),
          'time_interval'              => $self->o('time_interval'),
      },
      },
      { -logic_name => 'get_chr_dir',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
        -input_ids  => [],
        -parameters => {
          'input_dir' => $self->o('split_vcf_input_dir'),
          'inputcmd'  => 'ls #input_dir#',
        },
        -flow_into => { 
          '2->A' => {'init_spliceai' => {'input_chr_dir' => '#_0#'}},
          'A->1' => ['init_merge_files'],
        },
      },
      {   -logic_name => 'init_spliceai',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
          -input_ids  => [],
          -parameters => {
            'input_dir' => $self->o('split_vcf_input_dir') . '/' . '#input_chr_dir#',
            'inputcmd'  => 'find #input_dir# -type f -name "*.vcf" -printf "%f\n"',
          },
          -flow_into  => {
            2 => {'run_spliceai' => {'input_file' => '#_0#'}},
          }
      },
      { -logic_name => 'run_spliceai',
        -module => 'Bio::EnsEMBL::Variation::Pipeline::SpliceAI::RunSpliceAI',
        -hive_capacity => $self->o('pipeline_wide_analysis_capacity'),
        -analysis_capacity => $self->o('pipeline_wide_analysis_capacity'),
        -input_ids  => [],
        -rc_name => 'gpu',
        -parameters => {
          'main_dir'             => $self->o('main_dir'),
          'split_vcf_input_dir'  => $self->o('split_vcf_input_dir'),
          'split_vcf_output_dir' => $self->o('split_vcf_output_dir'),
          'fasta_file'           => $self->o('fasta_file'),
          'gene_annotation'      => $self->o('gene_annotation'),
          'masked_scores'        => $self->o('masked_scores'),
      },
        -max_retry_count => 3,
      },
      { -logic_name => 'init_merge_files',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
        -input_ids  => [],
        -parameters => {
          'input_dir' => $self->o('split_vcf_output_dir'),
          'inputcmd'  => 'ls #input_dir#',
        },
        -flow_into => { 
          2 => {'merge_files' => {'chr_dir' => '#_0#'}},
        },
      },
      { -logic_name => 'merge_files',
        -module => 'Bio::EnsEMBL::Variation::Pipeline::SpliceAI::MergeFiles',
        -input_ids  => [],
        -rc_name => '4Gb_job',
        -parameters => {
          'input_dir'        => $self->o('split_vcf_output_dir'),
          'output_dir'       => $self->o('output_dir'),
          'output_file_name' => $self->o('output_file_name'),
        },
        -flow_into => { 
          1 => ['finish_files'],
        },
      },
      { -logic_name => 'finish_files',
        -module => 'Bio::EnsEMBL::Variation::Pipeline::SpliceAI::FinishFiles',
        -input_ids  => [],
        -parameters => {
          'split_vcf_no_header_dir'   => $self->o('split_vcf_no_header_dir'),
          'split_vcf_input_dir'       => $self->o('split_vcf_input_dir'),
          'split_vcf_output_dir'      => $self->o('split_vcf_output_dir'),
        },
      }
  );
  return \@analyses;
}
1;
