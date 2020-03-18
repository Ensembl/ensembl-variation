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

package Bio::EnsEMBL::Variation::Pipeline::DBSNPImport::DBSNPImport_conf;

use strict;
use warnings;

use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
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
#
  return {
    %{ $self->SUPER::default_options() },   # inherit from parent
    # the general options that you should change to suit your environment
    hive_force_init         => 1,
    hive_use_param_stack    => 0,
    hive_use_triggers       => 0,
    hive_auto_rebalance_semaphores => 0,  # do not attempt to rebalance semaphores periodically by default
    hive_no_init            => 0, # setting it to 1 will skip pipeline_create_commands (useful for topping up)

    # the location of your checkout of the ensembl API
    hive_root_dir           => $ENV{'HOME'} . '/bin/ensembl-hive',
    ensembl_cvs_root_dir    => $ENV{'HOME'} . '/bin',

    debug                   => 0,

    pipeline_name           => 'dbsnp_import',
    species                 => 'homo_sapiens',
    assembly                => $self->o('assembly'),
    pipeline_dir            => $self->o('pipeline_dir'),
    data_dir                => $self->o('pipeline_dir') . '/split-src',
    rpt_dir                 => $self->o('pipeline_dir') . '/split-rpt',
    script_dir              => $self->o('ensembl_cvs_root_dir') . '/ensembl-variation/scripts/import/dbSNP_v2',
    ancestral_fasta_file    => $self->o('ancestral_fasta_file'),
    fasta_file              => $self->o('fasta_file'),
    registry_file           => $self->o('pipeline_dir') . '/' . 'ensembl.registry',

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
    %{$self->SUPER::pipeline_wide_parameters},          # here we inherit anything from the base class
    debug                => $self->o('debug'),
    
    pipeline_dir         => $self->o('pipeline_dir'),
    ensembl_registry     => $self->o('registry_file'),
    species              => $self->o('species'),
    assembly             => $self->o('assembly'),
    
    data_dir             => $self->o('data_dir'),
    rpt_dir              => $self->o('rpt_dir'),
    ancestral_fasta_file => $self->o('ancestral_fasta_file'),
    fasta_file           => $self->o('fasta_file'),
    script_dir           => $self->o('script_dir'),
  };
}

sub resource_classes {
    my ($self) = @_;
    return {
        %{$self->SUPER::resource_classes},  # inherit 'default' from the parent class
            'test_mem'    => { 'LSF' => '-R"select[mem>100] rusage[mem=100]" -M100'},
            'default_mem' => { 'LSF' => '-R"select[mem>1000] rusage[mem=1000]" -M1000'},
            'medium_mem'  => { 'LSF' => '-R"select[mem>4000] rusage[mem=4000]" -M4000'},
            'high_mem'    => { 'LSF' => '-R"select[mem>8000] rusage[mem=8000]" -M8000'},
    };
}

sub pipeline_analyses {
  my ($self) = @_;
  my @analyses;
  push @analyses, (
    {
      -logic_name => 'init_dbsnp_import',
      -module     => 'Bio::EnsEMBL::Variation::Pipeline::DBSNPImport::InitDBSNPImport',
      -input_ids  => [{}],
      -rc_name    => 'default_mem',
      -max_retry_count => 0,
      -flow_into  => {
        '2->A' => ['find_files'],
        'A->1' => ['report_dbsnp_import'],
      }
    },
    {
       -logic_name => 'find_files',
       -module => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
       -parameters => {
            'inputcmd' => 'ls #data_dir#/#sub_dir#',
            'column_names' => ['filename'],
        },
       -flow_into => {
          2 => {'load_dbsnp_file' => INPUT_PLUS()},
       },
    },
    {
      -logic_name        => 'load_dbsnp_file',
      -module            => 'Bio::EnsEMBL::Variation::Pipeline::DBSNPImport::LoadDBSNPFile',
      -rc_name           => 'default_mem',
      -max_retry_count   => 0,
      -analysis_capacity => 8,
    },
    {
      -logic_name => 'report_dbsnp_import',
      -module     => 'Bio::EnsEMBL::Variation::Pipeline::DBSNPImport::ReportDBSNPImport',
      -max_retry_count   => 0,
    },
  );
   return \@analyses;
}

1;
