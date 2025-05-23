=head1 LICENSE
Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2025] EMBL-European Bioinformatics Institute
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

package Bio::EnsEMBL::Variation::Pipeline::DumpHGVS::DumpHGVS_conf;

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
    ensembl_cvs_root_dir    => $ENV{'ENSEMBL_ROOT_DIR'} || $self->o('ensembl_cvs_root_dir'),
    hive_root_dir           => $self->o('ensembl_cvs_root_dir') . '/ensembl-hive',

    debug                   => 0,

    pipeline_name           => 'hgvs_dump',
    species                 => 'homo_sapiens',
    pipeline_dir            => $self->o('pipeline_dir'),
    hgvs_dir                => $self->o('pipeline_dir') . '/hgvs',
    registry_file           => $self->o('pipeline_dir') . '/' . 'ensembl.registry',

    region_size             => 5000000,
    region_overlap          =>    1000,
    bin_size                =>  500000,

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

sub pipeline_wide_parameters {
  my ($self) = @_;
  return {
    %{$self->SUPER::pipeline_wide_parameters},          # here we inherit anything from the base class
    debug            => $self->o('debug'),
    pipeline_dir     => $self->o('pipeline_dir'),
    ensembl_registry => $self->o('registry_file'),
    hgvs_dir         => $self->o('hgvs_dir'),
    species          => $self->o('species'),
  };
}

sub resource_classes {
    my ($self) = @_;
    return {
        %{$self->SUPER::resource_classes},  # inherit 'default' from the parent class
            'test_mem'    => { 'LSF' => '-q production -R"select[mem>100] rusage[mem=100]" -M100',
                               'SLURM' => '--partition=production --time=12:00:00 --mem=1G' },
            'default_mem' => { 'LSF' => '-q production -R"select[mem>2000] rusage[mem=2000]" -M2000',
                               'SLURM' => '--partition=production --time=12:00:00 --mem=4G' },
            'medium_mem'  => { 'LSF' => '-q production -R"select[mem>4000] rusage[mem=4000]" -M4000',
                               'SLURM' => '--partition=production --time=12:00:00 --mem=8G'},
            'high_mem'    => { 'LSF' => '-q production -R"select[mem>8000] rusage[mem=8000]" -M8000',
                               'SLURM' => '--partition=production --time=12:00:00 --mem=12G' },
    };
}

sub pipeline_analyses {
  my ($self) = @_;
  my @analyses;
  push @analyses, (
    {
      -logic_name => 'init_dump_hgvs',
      -module     => 'Bio::EnsEMBL::Variation::Pipeline::DumpHGVS::InitDumpHGVS',
      -parameters => {
        region_size => $self->o('region_size'),
        region_overlap => $self->o('region_overlap'),
      },
      -input_ids  => [{}],
      -rc_name    => 'default_mem',
      -max_retry_count => 0,
      -flow_into  => {
        '2->A' => ['dump_region_hgvs'],
        'A->1' => ['finish_dump_hgvs'],
      }
    },
    {
      -logic_name        => 'dump_region_hgvs',
      -module            => 'Bio::EnsEMBL::Variation::Pipeline::DumpHGVS::DumpRegionHGVS',
      -parameters     => {
          region_overlap  => $self->o('region_overlap'),
          bin_size => $self->o('bin_size')
      },
      -rc_name           => 'medium_mem',
      -max_retry_count   => 2,
      -analysis_capacity => 400,
    },
    {
      -logic_name => 'finish_dump_hgvs',
      -module     => 'Bio::EnsEMBL::Variation::Pipeline::DumpHGVS::FinishDumpHGVS',
      -rc_name    => 'default_mem',
      -max_retry_count   => 0,
    },
  );
   return \@analyses;
}

1;
