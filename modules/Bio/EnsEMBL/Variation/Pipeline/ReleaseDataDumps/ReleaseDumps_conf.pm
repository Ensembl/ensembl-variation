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
package Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::ReleaseDumps_conf;

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
        ensembl_release    => $self->o('ensembl_release'),

        # include or exclude the following species from the dumps, run for a division or all the species on the server
        species => [],
        antispecies => [],
        division    => [],
        run_all     => 1,

        pipeline_name      => $self->o('pipeline_name'),

        pipeline_dir       => $self->o('pipeline_dir'),
        
        ensembl_registry   => $self->o('ensembl_registry'),

        script_dir         => $self->o('ensembl_cvs_root_dir') . '/ensembl-variation/scripts',

        gvf_validator      => 'gvf_validator',

        vcf_validator      => 'vcf-validator',
        vcf_sort           => 'vcf-sort',

        so_file            => '/nfs/panda/ensembl/production/ensprod/obo_files/SO.obo',

        tmp_dir           => $self->o('pipeline_dir') . '/tmp_dir',
        gvf_readme => $self->o('ensembl_cvs_root_dir') . '/ensembl-variation/modules/Bio/EnsEMBL/Variation/Pipeline/ReleaseDataDumps/README_GVF',
        gvf_readme_human => $self->o('ensembl_cvs_root_dir') . '/ensembl-variation/modules/Bio/EnsEMBL/Variation/Pipeline/ReleaseDataDumps/README_GVF_human',

        vcf_readme => $self->o('ensembl_cvs_root_dir') . '/ensembl-variation/modules/Bio/EnsEMBL/Variation/Pipeline/ReleaseDataDumps/README_VCF',
        vcf_readme_human => $self->o('ensembl_cvs_root_dir') . '/ensembl-variation/modules/Bio/EnsEMBL/Variation/Pipeline/ReleaseDataDumps/README_VCF_human',

        global_vf_count_in_species => 5_000_000, # if number of vf in a species exceeds this we need to split up dumps
        max_vf_load => 2_000_000, # group slices together until the vf count exceeds max_vf_load
        vf_per_slice => 2_000_000, # if number of vf exceeds this we split the slice and dump for each split slice
        max_split_slice_length => 5e6, # 1e7

        data_dir => '/nfs/production/panda/ensembl/variation/data/',

        ancestral_alleles_file_dir => {'homo_sapiens' => {
                                          'GRCh37' => $self->o('data_dir') . 'ancestral_alleles/GRCh37',
                                          'GRCh38' => $self->o('data_dir') . 'ancestral_alleles/' . $self->o('ensembl_release'),
                                        },
                                      },
        fasta_file  => {'homo_sapiens' => {
                              'GRCh37' => $self->o('data_dir') . 'Homo_sapiens.GRCh37.75.dna.primary_assembly.fa',
                              'GRCh38' => $self->o('data_dir') . 'Homo_sapiens.GRCh38.dna.toplevel.fa.gz',
                            },
                          },

        debug => 0,

        # init_pipeline.pl will create the hive database on this machine, naming it
        # <username>_<pipeline_name>, and will drop any existing database with this
        # name
       
        pipeline_wide_analysis_capacity => 50,        
        pipeline_wide_hive_capacity => 50,        
 
        only_finish_dumps => 0,
        human_population_dumps => 0,

        pipeline_db => {
            -host   => $self->o('hive_db_host'),
            -port   => $self->o('hive_db_port'),
            -user   => $self->o('hive_db_user'),
            -pass   => $self->o('hive_db_password'),            
            -dbname => $ENV{'USER'} . '_' . $self->o('pipeline_name'),
            -driver => 'mysql',
        },
    };
}

sub pipeline_wide_parameters {
    my ($self) = @_;
    return {
        %{$self->SUPER::pipeline_wide_parameters}, # here we inherit anything from the base class
        release          => $self->o('ensembl_release'),
        pipeline_dir     => $self->o('pipeline_dir'),
        script_dir       => $self->o('script_dir'), 
        so_file          => $self->o('so_file'),    
        gvf_validator    => $self->o('gvf_validator'),
        vcf_validator    => $self->o('vcf_validator'),
        vcf_sort         => $self->o('vcf_sort'),
        tmp_dir          => $self->o('tmp_dir'),
        gvf_readme       => $self->o('gvf_readme'), 
        vcf_readme       => $self->o('vcf_readme'),
        gvf_readme_human => $self->o('gvf_readme_human'), 
        vcf_readme_human => $self->o('vcf_readme_human'),
        pipeline_wide_analysis_capacity => $self->o('pipeline_wide_analysis_capacity'),        
        debug => $self->o('debug'),
        global_vf_count_in_species => $self->o('global_vf_count_in_species'),
        max_vf_load => $self->o('max_vf_load'),
        vf_per_slice => $self->o('vf_per_slice'),
        max_split_slice_length => $self->o('max_split_slice_length'),
        division => $self->o('division'),
        ensembl_registry => $self->o('ensembl_registry'),
    };
}

# Override the default method, to force an automatic loading of the registry in all workers
sub beekeeper_extra_cmdline_options {
  my ($self) = @_;
  return
      ' -reg_conf ' . $self->o('ensembl_registry'),
  ;
}

sub resource_classes {
    my ($self) = @_;
    return {
        %{$self->SUPER::resource_classes},
        'default' => { 'LSF' => '-q production-rh74 -R"select[mem>1500] rusage[mem=1500]" -M1500'},
        'medium'  => { 'LSF' => '-q production-rh74 -R"select[mem>4500] rusage[mem=4500]" -M4500'},
        'high'    => { 'LSF' => '-q production-rh74 -R"select[mem>8500] rusage[mem=8500]" -M8500'},
    };
}
sub pipeline_analyses {
  my ($self) = @_;
  my @analyses;
  push @analyses, (
      {   -logic_name => 'pre_run_checks',
          -module => 'Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::PreRunChecks',
          -max_retry_count => 0,
          -input_ids     => [{}],
          -flow_into => {
              1 => ['species_factory'],
          },
          -parameters => {
              'file_type' => 'gvf',
          },
      },
      {
         -logic_name    => 'species_factory',
         -module         => 'Bio::EnsEMBL::Production::Pipeline::Common::SpeciesFactory',
         -parameters    => {
           species     => $self->o('species'),
           antispecies => $self->o('antispecies'),
           division    => $self->o('division'),
           run_all     => $self->o('run_all'),
          },
         -rc_name       => 'default',
         -hive_capacity => 1,
         -max_retry_count => 0,
         -flow_into     => {
             4 => WHEN(
              '#only_finish_dumps#' => 'finish_dumps',
              ELSE 'generate_config',
            )
          },
      },
      {   -logic_name => 'generate_config',
          -module => 'Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::Config',
          -max_retry_count => 0,
          -rc_name => 'default',
          -flow_into     => {
          '2->A' => ['init_dump'],
          'A->1' => ['init_join_split_slice']
        }
      },
      {   -logic_name => 'init_dump',
          -module => 'Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::InitSubmitJob',
          -max_retry_count => 0,
          -analysis_capacity => $self->o('pipeline_wide_analysis_capacity'),
          -flow_into => {
              2 => ['submit_job_gvf_dumps'],
          },
          -parameters => {
              'file_type' => 'gvf',
              'job_type'  => 'dump',
          },
          -rc_name => 'default',
      },
      {   -logic_name => 'submit_job_gvf_dumps',
          -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -analysis_capacity  => $self->o('pipeline_wide_analysis_capacity'),
          -hive_capacity => $self->o('pipeline_wide_analysis_capacity'),
          -max_retry_count => 1,
          -rc_name => 'default',
          -flow_into      => {
            -1 => ['submit_job_gvf_dumps_mediummem'],
          }
      },
      {   -logic_name => 'submit_job_gvf_dumps_mediummem',
          -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -analysis_capacity  => $self->o('pipeline_wide_analysis_capacity'),
          -hive_capacity => $self->o('pipeline_wide_analysis_capacity'),
          -max_retry_count => 1,
          -rc_name => 'medium',
          -flow_into      => {
            -1 => ['submit_job_gvf_dumps_highmem'],
          }
      },
      {   -logic_name => 'submit_job_gvf_dumps_highmem',
          -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -analysis_capacity  => $self->o('pipeline_wide_analysis_capacity'),
          -hive_capacity => $self->o('pipeline_wide_analysis_capacity'),
          -max_retry_count => 0,
          -rc_name => 'high',
      },
# join split slice
      { -logic_name => 'init_join_split_slice',
        -module => 'Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::InitJoinDump',
        -analysis_capacity => $self->o('pipeline_wide_analysis_capacity'),
        -parameters => {
            'mode' => 'join_split_slice',
            'file_type' => 'gvf',
          },
        -flow_into => {
          '2->A' => ['join_split_slice'],
          'A->1' => ['init_validate_gvf'],
        },
      },
      { -logic_name => 'join_split_slice',
        -module => 'Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::JoinDump',
        -analysis_capacity => $self->o('pipeline_wide_analysis_capacity'),
      },
# validate 
      { -logic_name => 'init_validate_gvf',
        -module => 'Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::InitValidation',
        -analysis_capacity => $self->o('pipeline_wide_analysis_capacity'),
        -flow_into => {
          '2->A' => ['validate_gvf'],
          'A->1' => ['summary_validate_gvf'],
        },
        -parameters => {
          'file_type' => 'gvf',
        },
      },
      { -logic_name => 'validate_gvf',
        -module => 'Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::Validate',
        -parameters => {
          'file_type'     => 'gvf',
          'so_file'       => $self->o('so_file'),
          'gvf_validator' => $self->o('gvf_validator'),
        },
      },
      { -logic_name => 'summary_validate_gvf', # die if Error?
        -module => 'Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::SummariseValidation',
        -parameters => {
          'file_type' => 'gvf',
        },
        -flow_into => {
          1 => ['cleanup_gvf_dumps'],
        },
      },

# cleanup gvf dumps
      { -logic_name        => 'cleanup_gvf_dumps',
        -analysis_capacity => $self->o('pipeline_wide_analysis_capacity'),
        -module            => 'Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::FileUtils',
        -parameters        => {
          'mode'  => 'post_gvf_dump_cleanup',
        },
        -flow_into => {
            '2->A' => ['init_parse'],
            'A->1' => ['summary_validate_vcf']
        },
      },

# gvf2vcf
      {   -logic_name => 'init_parse',
          -module => 'Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::InitSubmitJob',
          -flow_into => {
              2 => ['submit_job_gvf2vcf'],
          },
          -parameters => {
              'file_type' => 'vcf',
              'job_type'  => 'parse',
              'ancestral_alleles_file_dir' => $self->o('ancestral_alleles_file_dir'),
              'fasta_file' => $self->o('fasta_file'),

          },
      },
      {   -logic_name => 'submit_job_gvf2vcf',
          -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -analysis_capacity  => $self->o('pipeline_wide_analysis_capacity'),
          -max_retry_count => 0,
          -rc_name => 'default',
          -flow_into => ['validate_vcf'],
      },
      {   -logic_name => 'validate_vcf',
          -module => 'Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::Validate',
          -parameters => {
              'file_type' => 'vcf',
          },
      },
      { -logic_name => 'summary_validate_vcf',
        -module => 'Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::SummariseValidation',
          -parameters => {
              'file_type' => 'vcf',
          },
          -flow_into => {
              1 => ['cleanup_gvf2vcf'],
          },
      },

      {   -logic_name => 'cleanup_gvf2vcf',
          -module     => 'Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::FileUtils',
          -parameters => {
            'mode'     => 'post_gvf2vcf_cleanup',
          },
          -flow_into => {
              '2->A' => ['init_join_dumps'],
              'A->1' => ['cleanup_dumps']
           },
      },
     {   -logic_name => 'init_join_dumps',
          -module => 'Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::InitJoinDump',
          -flow_into => {
              2 => ['join_dumps'],
          },
      },
      {   -logic_name => 'join_dumps',
          -module => 'Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::JoinDump',
      },
      {   -logic_name => 'cleanup_dumps',
          -module => 'Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::CleanUp',
          -parameters => {
            'mode' => 'post_join_dumps',
          },
          -flow_into => {
              1 => ['finish_dumps'],
          },
      },
      {   -logic_name => 'finish_dumps',
          -module => 'Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::Finish',
          -parameters => {
              'gvf_readme' => $self->o('gvf_readme'),
              'vcf_readme' => $self->o('vcf_readme'),
              'gvf_readme_human' => $self->o('gvf_readme_human'),
              'vcf_readme_human' => $self->o('vcf_readme_human'),
          }
      },
  );
  return \@analyses;
}
1;
