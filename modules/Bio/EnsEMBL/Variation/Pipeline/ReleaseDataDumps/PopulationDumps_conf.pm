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
package Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::PopulationDumps_conf;

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
        hive_use_param_stack => 0,
        hive_use_triggers    => 0,
        hive_root_dir        => $ENV{'HOME'} . '/bin/ensembl-hive', 
        ensembl_cvs_root_dir => $ENV{'HOME'} . '/bin',
        hive_no_init         => 0,


        ensembl_release    => $self->o('ensembl_release'),
        pipeline_name      => $self->o('pipeline_name'),

        pipeline_dir       => $self->o('pipeline_dir'),
        
        registry_file      => $self->o('registry_file'),

        script_dir         => $self->o('ensembl_cvs_root_dir') . '/ensembl-variation/scripts',

        gvf_validator      => 'gvf_validator',
        so_file            => '/nfs/panda/ensembl/production/ensprod/obo_files/SO.obo',

        tmp_dir           => $self->o('tmp_dir'),
        gvf_readme => $self->o('ensembl_cvs_root_dir') . '/ensembl-variation/modules/Bio/EnsEMBL/Variation/Pipeline/ReleaseDataDumps/README_GVF',
        vcf_readme => $self->o('ensembl_cvs_root_dir') . '/ensembl-variation/modules/Bio/EnsEMBL/Variation/Pipeline/ReleaseDataDumps/README_VCF',

        # init_pipeline.pl will create the hive database on this machine, naming it
        # <username>_<pipeline_name>, and will drop any existing database with this
        # name
        
        finish_dumps => 0,
        human_population_dumps => 0,

        hive_db_host    => 'ens-variation',
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
    };
}

sub resource_classes {
    my ($self) = @_;
    return {
        %{$self->SUPER::resource_classes},
        'default' => { 'LSF' => '-R"select[mem>5500] rusage[mem=5500]" -M5500'},
        'urgent'  => { 'LSF' => '-q yesterday -R"select[mem>2000] rusage[mem=2000]" -M2000'},
        'highmem' => { 'LSF' => '-R"select[mem>15000] rusage[mem=15000]" -M15000'}, # this is Sanger LSF speak for "give me 15GB of memory"
        'long'    => { 'LSF' => '-q long -R"select[mem>2000] rusage[mem=2000]" -M2000'},
    };
}

sub pipeline_analyses {
  my ($self) = @_;
  my @analyses;

  push @analyses, (
    { -logic_name => 'init_population_gvf',
      -module => 'Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::InitDumpPopulation',
      -input_ids => [{},],
      -parameters => {
        'file_type' => 'gvf',
        'prefetched_frequencies' => $self->o('prefetched_frequencies'),
      },
      -flow_into => {
        '2->A' => ['dump_population_gvf'],
        'A->1' => ['finish_dump_population_gvf']
      },
    },
    { -logic_name => 'dump_population_gvf',
      -module => 'Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::DumpPopulation',
      -parameters => {
        'file_type' => 'gvf',
        'prefetched_frequencies' => $self->o('prefetched_frequencies'),
      },
    },
    { -logic_name => 'finish_dump_population_gvf',
      -module => 'Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::FinishDumpPopulation',
      -parameters => {
        'file_type' => 'gvf',
      },
      -flow_into => {
        1 => ['init_population_vcf'],
      },
    },
    { -logic_name => 'init_population_vcf',
      -module => 'Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::InitDumpPopulation',
      -parameters => {
        'file_type' => 'vcf',
      },
      -flow_into => {
        '2->A' => ['dump_population_vcf'],
        'A->1' => ['finish_dump_population_vcf']
      },
    },
    { -logic_name => 'dump_population_vcf',
      -module => 'Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::SubmitJob',
    },
    { -logic_name => 'finish_dump_population_vcf',
      -module => 'Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::FinishDumpPopulation',
      -parameters => {
        'file_type' => 'vcf',
      },
    },
  );
  return \@analyses;
}

1;

