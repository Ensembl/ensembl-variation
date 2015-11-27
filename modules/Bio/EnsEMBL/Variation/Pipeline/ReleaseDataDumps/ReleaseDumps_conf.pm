=head1 LICENSE

Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');

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
        hive_auto_rebalance_semaphores => 1,
        hive_force_init      => 1,
        hive_use_param_stack => 1,
        hive_use_triggers    => 0,
        hive_root_dir        => $ENV{'HOME'} . '/DEV/ensembl-hive', 
        ensembl_cvs_root_dir => $ENV{'HOME'} . '/DEV',
        hive_no_init         => 0,


        ensembl_release    => $self->o('ensembl_release'),
        pipeline_name      => $self->o('pipeline_name'),

        pipeline_dir       => $self->o('pipeline_dir'),
        
        registry_file      => $self->o('registry_file'),

        script_dir         => $self->o('ensembl_cvs_root_dir') . '/ensembl-variation/scripts',

        gvf_validator      => '/nfs/users/nfs_a/at7/tools/gvf_validator',
        so_file            => '/nfs/users/nfs_a/at7/obo/e' . $self->o('ensembl_release') . '/so.obo',

        tmp_dir           => $self->o('tmp_dir'),
        gvf_readme => $self->o('ensembl_cvs_root_dir') . '/ensembl-variation/modules/Bio/EnsEMBL/Variation/Pipeline/ReleaseDataDumps/README_GVF',
        vcf_readme => $self->o('ensembl_cvs_root_dir') . '/ensembl-variation/modules/Bio/EnsEMBL/Variation/Pipeline/ReleaseDataDumps/README_VCF',

        # init_pipeline.pl will create the hive database on this machine, naming it
        # <username>_<pipeline_name>, and will drop any existing database with this
        # name
       
        pipeline_wide_analysis_capacity => 50,        
 
        only_finish_dumps => 0,
        human_population_dumps => 0,

        hive_db_host    => 'ens-variation3',
        hive_db_port    => 3306,
        hive_db_user    => 'ensadmin',

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
        registry_file    => $self->o('registry_file'),
        output_path      => $self->o('pipeline_dir'),
        pipeline_dir     => $self->o('pipeline_dir'),
        script_dir       => $self->o('script_dir'), 
        so_file          => $self->o('so_file'),    
        gvf_validator    => $self->o('gvf_validator'),
        tmp_dir          => $self->o('tmp_dir'),
        gvf_readme       => $self->o('gvf_readme'), 
        vcf_readme       => $self->o('vcf_readme'),
        pipeline_wide_analysis_capacity => $self->o('pipeline_wide_analysis_capacity'),        

        config_file => '/lustre/scratch110/ensembl/at7/release_83/dumps/data_dumps_config.json',
    };
}

sub resource_classes {
    my ($self) = @_;
    return {
        %{$self->SUPER::resource_classes},
        'default' => { 'LSF' => '-q normal -R"select[mem>3500] rusage[mem=3500]" -M3500'},
        'urgent'  => { 'LSF' => '-q yesterday -R"select[mem>2000] rusage[mem=2000]" -M2000'},
        'highmem' => { 'LSF' => '-q long -R"select[mem>15000] rusage[mem=15000]" -M15000'}, # this is Sanger LSF speak for "give me 15GB of memory"
        'long'    => { 'LSF' => '-q long -R"select[mem>2000] rusage[mem=2000]" -M2000'},
    };
}
sub pipeline_analyses {
  my ($self) = @_;
  my @analyses;
  if ($self->o('only_finish_dumps')) {
    push @analyses, (
      { -logic_name => 'finish_dumps',
        -module => 'Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::Finish',
        -input_ids => [{},],
      },
    );
  } else {
    push @analyses, (
        {   -logic_name => 'pre_run_checks_gvf_dumps',
            -module => 'Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::PreRunChecks',
            -input_ids => [{},],
            -max_retry_count => 1,
            -flow_into => {
                1 => ['export_config_file'],
            },
            -parameters => {
                'file_type' => 'gvf',   
            },
        },
        {   -logic_name => 'export_config_file',
            -module => 'Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::Config',
            -max_retry_count => 1,
            -rc_name => 'default',
             -flow_into => {
                1 => ['species_factory_gvf_dumps'],
            },
        },
        {   -logic_name => 'species_factory_gvf_dumps',
            -module => 'Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::SpeciesFactory', 
            -analysis_capacity  => $self->o('pipeline_wide_analysis_capacity'),
            -flow_into => {
                '2->A' => ['init_dump'],
                'A->1' => ['species_factory_join_split_slice'],
             },
        },
        {   -logic_name => 'init_dump',
            -module => 'Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::InitSubmitJob',
            -max_retry_count => 1,
            -analysis_capacity => $self->o('pipeline_wide_analysis_capacity'),
            -flow_into => {
                1 => ['submit_job_gvf_dumps'],
            },
            -parameters => {
                'file_type' => 'gvf',
                'job_type'  => 'dump',
            },
            -rc_name => 'default',
        },
        {   -logic_name => 'submit_job_gvf_dumps',
            -module => 'Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::SubmitJob',
            -analysis_capacity  => $self->o('pipeline_wide_analysis_capacity'),
            -max_retry_count => 1,
            -rc_name => 'default',
        },


# join split slice
        { -logic_name => 'species_factory_join_split_slice',
          -module => 'Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::SpeciesFactory',
          -flow_into => {
            '2->A' => ['join_split_slice'],
            'A->1' => ['species_factory_validate_gvf'],
          },
        },
        { -logic_name => 'join_split_slice',
          -module => 'Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::JoinDump',
          -analysis_capacity => $self->o('pipeline_wide_analysis_capacity'),
          -parameters => {
              'mode' => 'join_slice_split',
              'file_type' => 'gvf',
            },
        },


# validate 
        { -logic_name => 'species_factory_validate_gvf',
          -module => 'Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::SpeciesFactory',
          -flow_into => {
            '2->A' => ['init_validate_gvf'],
            'A->1' => ['summary_validate_gvf']
          },
          -parameters => {
            'file_type' => 'gvf',
          },
        },
        { -logic_name => 'init_validate_gvf',
          -module => 'Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::InitValidation',
          -analysis_capacity => $self->o('pipeline_wide_analysis_capacity'),
          -flow_into => {
              1 => ['validate_gvf'],
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
            1 => ['species_factory_cleanup_gvf_dumps'],
          },
        },


# cleanup gvf dumps
        {   -logic_name => 'species_factory_cleanup_gvf_dumps',
            -module => 'Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::SpeciesFactory',
            -flow_into => {
             '2->A' => ['cleanup_gvf_dumps'],
             'A->1' => ['pre_run_checks_gvf2vcf' ],
            },
        },
        { -logic_name        => 'cleanup_gvf_dumps',
          -analysis_capacity => $self->o('pipeline_wide_analysis_capacity'),
          -module            => 'Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::FileUtils',
          -parameters        => {
            'mode'  => 'post_gvf_dump_cleanup',
          },
        },

# gvf2vcf
        {   -logic_name => 'pre_run_checks_gvf2vcf',
            -module => 'Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::PreRunChecks',
            -max_retry_count => 1,
            -flow_into => {
                1 => ['species_factory_gvf2vcf'],
            },
            -parameters => {
                'file_type' => 'vcf',
            },
        },
        {   -logic_name => 'species_factory_gvf2vcf',
            -module => 'Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::SpeciesFactory',
            -analysis_capacity => 3,
            -flow_into => {
                '2->A' => ['init_parse'],
                'A->1' => ['summary_validate_vcf']
             },
            -parameters => {
                'file_type' => 'vcf',
            },
        },
        {   -logic_name => 'init_parse',
            -module => 'Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::InitSubmitJob',
            -flow_into => {
                1 => ['submit_job_gvf2vcf'],
            },
            -parameters => {
                'file_type' => 'vcf',
                'job_type'  => 'parse',
            },
        },
        {   -logic_name => 'submit_job_gvf2vcf',
            -module => 'Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::SubmitJob',
            -analysis_capacity => $self->o('pipeline_wide_analysis_capacity'),
            -flow_into => {
                1 => ['validate_vcf'],
            },
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
              1 => ['species_factory_cleanup_gvf2vcf'],
            },
        },

        {   -logic_name => 'species_factory_cleanup_gvf2vcf',
            -module => 'Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::SpeciesFactory',
            -analysis_capacity => 3,
            -flow_into => {
                '2->A' => ['cleanup_gvf2vcf'],
                'A->1' => ['species_factory_join_dumps']
             },
        },
        {   -logic_name => 'cleanup_gvf2vcf',
            -module     => 'Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::FileUtils',
            -parameters => {
              'mode'     => 'post_gvf2vcf_cleanup',
            },
        },
        {   -logic_name => 'species_factory_join_dumps',
            -module => 'Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::SpeciesFactory', 
            -analysis_capacity  => $self->o('pipeline_wide_analysis_capacity'),
            -flow_into => {
                '2->A' => ['init_join_dumps'],
                'A->1' => ['species_factory_cleanup_dumps']
             },
        },
       {   -logic_name => 'init_join_dumps',
            -module => 'Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::InitJoinDump',
            -flow_into => {
                1 => ['join_dumps'],
            },
        },
        {   -logic_name => 'join_dumps',
            -module => 'Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::JoinDump',
        },
        {   -logic_name => 'species_factory_cleanup_dumps',
            -module => 'Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::SpeciesFactory', 
            -analysis_capacity  => $self->o('pipeline_wide_analysis_capacity'),
            -flow_into => {
                '2->A' => ['cleanup_dumps'],
                'A->1' => ['finish_dumps']
             },
        },
        {   -logic_name => 'cleanup_dumps',
            -module => 'Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::CleanUp',
            -parameters => {
              'mode' => 'post_join_dumps',
            },
        },
        {   -logic_name => 'finish_dumps',
            -module => 'Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::Finish',
        },
    );
  }
  return \@analyses;
}
1;
