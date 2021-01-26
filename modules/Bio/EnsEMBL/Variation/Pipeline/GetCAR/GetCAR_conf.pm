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

package Bio::EnsEMBL::Variation::Pipeline::GetCAR::GetCAR_conf;

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

    pipeline_name           => 'car_lu',
    species                 => 'homo_sapiens',
    pipeline_dir            => $self->o('pipeline_dir'),
    hgvs_dir                => $self->o('pipeline_dir') . '/hgvs',
    car_lu_dir              => $self->o('pipeline_dir') . '/car_lu',

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
    debug            => $self->o('debug'),
    pipeline_dir     => $self->o('pipeline_dir'),
    hgvs_dir         => $self->o('hgvs_dir'),
    car_lu_dir       => $self->o('car_lu_dir'),
    species          => $self->o('species'),
  };
}

sub resource_classes {
    my ($self) = @_;
    return {
        %{$self->SUPER::resource_classes},  # inherit 'default' from the parent class
            'test_mem'    => { 'LSF' => '-R"select[mem>100] rusage[mem=100]" -M100'},
            'default_mem' => { 'LSF' => '-R"select[mem>1000] rusage[mem=1000]" -M1000'},
            'medium_mem'  => { 'LSF' => '-R"select[mem>2000] rusage[mem=2000]" -M2000'},
            'high_mem'    => { 'LSF' => '-R"select[mem>4000] rusage[mem=4000]" -M4000'},
    };
}

sub pipeline_analyses {
  my ($self) = @_;
  my @analyses;
  push @analyses, (
    {
      -logic_name => 'init_get_car',
      -module     => 'Bio::EnsEMBL::Variation::Pipeline::GetCAR::InitGetCAR',
      -parameters => {
        hgvs_dir => $self->o('hgvs_dir'),
        car_lu_dir => $self->o('car_lu_dir'),
      } ,
      -input_ids => [{}],
      -rc_name    => 'default_mem',
      -max_retry_count => 0,
      -flow_into  => {
        '2->A' => ['get_car_seqname'],
        'A->1' => ['finish_get_car'],
      }
    },
    {
      -logic_name        => 'get_car_seqname',
      -module            => 'Bio::EnsEMBL::Variation::Pipeline::GetCAR::GetCARSeqname',
      -parameters     => {
          hgvs_dir => $self->o('hgvs_dir'),
          car_lu_dir => $self->o('car_lu_dir'),
      },
      -rc_name           => 'default_mem',
      -max_retry_count   => 0,
      -analysis_capacity => 1,
    },
    {
      -logic_name => 'finish_get_car',
      -module     => 'Bio::EnsEMBL::Variation::Pipeline::GetCAR::FinishGetCAR',
      -max_retry_count   => 0,
    },
  );
   return \@analyses;
}

1;
