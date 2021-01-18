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
 developers list at <https://lists.ensembl.org/mailman/listinfo/dev>.
 Questions may also be sent to the Ensembl help desk at
 <https://www.ensembl.org/Help/Contact>.
=cut
package Bio::EnsEMBL::Variation::Pipeline::RemappingVCF::RemappingVCF_conf;

use strict;
use warnings;

use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
use base ('Bio::EnsEMBL::Hive::PipeConfig::EnsemblGeneric_conf');

sub default_options {
    my ($self) = @_;
    return {
        %{ $self->SUPER::default_options()
        },

        hive_auto_rebalance_semaphores   => 1,
        hive_force_init                  => 1,
        hive_use_param_stack             => 1,
        population                       => $self->o('population'),
        pipeline_dir                     => $self->o('pipeline_dir'),
        registry_file_newasm             => $self->o('pipeline_dir') . 'ensembl.registry.newasm',
        registry_file_oldasm             => $self->o('pipeline_dir') . 'ensembl.registry.oldasm',
        registry_file_oldasm_same_server => $self->o('pipeline_dir') . 'ensembl.registry.oldasm.same_server',
        vcf_file                         => $self->o('vcf_file'),
        pipeline_name                    => $self->o('pipeline_name'), 
        load_from_vcf                    => 1 ,

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
        pipeline_dir  => $self->o('pipeline_dir'),
        registry_file_newasm => $self->o('registry_file_newasm'),
        registry_file_oldasm => $self->o('registry_file_oldasm'),
        registry_file_oldasm_same_server => $self->o('registry_file_oldasm_same_server'),
        species       => $self->o('species'),
        vcf_file      => $self->o('vcf_file'),
        population    => $self->o('population'),
    };
}

sub resource_classes {
    my ($self) = @_;
    return {
        %{$self->SUPER::resource_classes},
        'default' => { 'LSF' => '-q production-rh7 -R"select[mem>5500] rusage[mem=5500]" -M5500'},
    };
}

sub pipeline_analyses {
  my ($self) = @_;
  my @analyses;

  if ($self->o('load_from_vcf')) {
    push @analyses, (
      {
        -logic_name => 'pre_run_checks',
        -module     => 'Bio::EnsEMBL::IntVar::Pipeline::NextGen::RemappingVCF::PreRunChecks',
        -input_ids  => [{},],
        -max_retry_count => 0,
        -flow_into  => {
          1 => ['load_from_vcf']
        },
      },
      {
        -logic_name => 'load_from_vcf',
        -module     => 'Bio::EnsEMBL::IntVar::Pipeline::NextGen::RemappingVCF::LoadFromVCF',
        -rc_name => 'default',
        -max_retry_count => 0,
        -parameters => {
          'dump_data_from_VCF' => 1,
          'load_data' => 1,
          'update_mappings' => 1,
        },
        -flow_into  => {
          1 => ['init_remapping']
        },
      },
      {
         -logic_name => 'init_remapping',
         -module     => 'Bio::EnsEMBL::IntVar::Pipeline::NextGen::RemappingVCF::InitRemapping',
         -rc_name    => 'default',
        -max_retry_count => 0,
         -flow_into  => {
          '2->A' => ['remapping'],
          'A->1' => ['join_remapping'],
        },
      },
    );
  } else {
    push @analyses, (
      {
        -logic_name => 'pre_run_checks',
        -module     => 'Bio::EnsEMBL::IntVar::Pipeline::NextGen::RemappingVCF::PreRunChecks',
        -input_ids  => [{},],
        -max_retry_count => 0,
        -flow_into  => {
          1 => ['init_remapping']
        },
      },
      {
        -logic_name => 'init_remapping',
        -module     => 'Bio::EnsEMBL::IntVar::Pipeline::NextGen::RemappingVCF::InitRemapping',
        -input_ids  => [{},],
        -rc_name    => 'default',
        -max_retry_count => 0,
        -flow_into  => {
          '2->A' => ['remapping'],
          'A->1' => ['join_remapping'],
        },
      },
    );
  }
  push @analyses, (
      {   -logic_name => 'remapping',
          -module     => 'Bio::EnsEMBL::IntVar::Pipeline::NextGen::RemappingVCF::Remapping',
          -rc_name    => 'default',
          -hive_capacity  => 14,
          -max_retry_count => 0,
      },
      {   -logic_name => 'join_remapping',
          -module     => 'Bio::EnsEMBL::IntVar::Pipeline::NextGen::RemappingVCF::JoinRemapping',
          -rc_name    => 'default',
      },
  );
  return \@analyses;
}
1;
