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

package Bio::EnsEMBL::Variation::Pipeline::RegulationEffect_conf;

use strict;
use warnings;
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
use base ('Bio::EnsEMBL::Hive::PipeConfig::EnsemblGeneric_conf');

sub default_options {
    my ($self) = @_;
    return {
        %{$self->SUPER::default_options},
        # general pipeline options that you should change to suit your environment
        hive_default_max_retry_count   => 0,
        hive_force_init                => 1,
        hive_use_param_stack           => 0,
        hive_use_triggers              => 0, 
        hive_auto_rebalance_semaphores => 0,
        hive_no_init                   => 0,
        # the location of your checkout of the ensembl API (the hive looks for SQL files here)
        ensembl_cvs_root_dir           => $self->o('ensembl_cvs_root_dir'),
        hive_root_dir                  => $self->o('ensembl_cvs_root_dir') . '/ensembl-hive',

        pipeline_name                  => $self->o('pipeline_name'),
        species                        => $self->o('species'),
        pipeline_dir                   => $self->o('pipeline_dir'),
        registry_file                  => $self->o('pipeline_dir') . '/ensembl.registry',

        debug                           => 0,
        split_slice                     => 1,
        split_slice_length              => 5e6,
        only_update_vf                  => 0,

        # if set to 1 this option tells the transcript_effect analysis to disambiguate
        # ambiguity codes in single nucleotide alleles, so e.g. an allele string like
        # 'T/M' will be treated as if it were 'T/A/C' (this was a request from ensembl
        # genomes and we don't use it by default in the ensembl variation pipeline)
        disambiguate_single_nucleotide_alleles => 1,

        pipeline_db => {
            -host   => $self->o('hive_db_host'),
            -port   => $self->o('hive_db_port'),
            -user   => $self->o('hive_db_user'),
            -pass   => $self->o('hive_db_password'),
            -dbname => $ENV{'USER'} . '_ehive_' . $self->o('pipeline_name') . '_' . $self->o('ensembl_release') . '_' . $self->o('assembly') . '_' . $self->o('species'),
            -driver => 'mysql',
        },
    };
}


sub resource_classes {
    my ($self) = @_;
    # configuration for the various resource options used in the pipeline
    # EBI farm users should either change these here, or override them on the
    # command line to suit the EBI farm. The names of each option hopefully
    # reflect their usage, but you may want to change the details (memory
    # requirements, queue parameters etc.) to suit your own data

    return {
      'default' => { 'LSF' => '-qproduction -R"select[mem>2000] rusage[mem=2000]" -M2000',
                     'SLURM' => "--partition=production --time=4:00:00 --mem=2G" },
      'highmem' => { 'LSF' => '-qproduction -R"select[mem>15000] rusage[mem=15000]" -M15000',
                     'SLURM' => "--partition=production --time=28:00:00 --mem=15G"},
      'long'    => { 'LSF' => '-qproduction -R"select[mem>2000] rusage[mem=2000]" -M2000',
                     'SLURM' => "--partition=production --time=72:00:00 --mem=4G"},
    };
}

sub pipeline_wide_parameters {
    my ($self) = @_;
    return {
        %{$self->SUPER::pipeline_wide_parameters},
        disambiguate_single_nucleotide_alleles => $self->o('disambiguate_single_nucleotide_alleles'),
        ensembl_registry                       => $self->o('registry_file'),
        only_update_vf                         => $self->o('only_update_vf'),
        species                                => $self->o('species'),
        debug                                  => $self->o('debug'),
        split_slice                            => $self->o('split_slice'),
        split_slice_length                     => $self->o('split_slice_length'),
    };
}

sub pipeline_analyses {
    my ($self) = @_;
   
    my @analyses;
    if ($self->o('only_update_vf')) {
        push @analyses, (
            {   -logic_name => 'finish_regulation_effect',
                -module => 'Bio::EnsEMBL::Variation::Pipeline::FinishRegulationEffect',
                -input_ids => [{},],
                -hive_capacity => 1,
                -parameters => {
                  'pipeline_dir' => $self->o('pipeline_dir'),
                },
            },
        );
    } else {
        push @analyses, (
            {   -logic_name => 'init_regulation_effect',
                -module => 'Bio::EnsEMBL::Variation::Pipeline::InitRegulationEffect',
                -hive_capacity => 1,
                -rc_name => 'default',
                -flow_into => { 1 => ['init_slice'] },
                -rc_name => 'default',
                -input_ids => [{},],
            }, 
            {   -logic_name => 'init_slice',
                -module => 'Bio::EnsEMBL::Variation::Pipeline::SliceFactory',
                -hive_capacity => 1,
                -rc_name => 'default',
                -max_retry_count => 0,
                -flow_into => {
                    '2->A' => ['regulation_effect'],
                    'A->1' => ['finish_regulation_effect'],
                },
            },
            {   -logic_name => 'regulation_effect',
                -module => 'Bio::EnsEMBL::Variation::Pipeline::RegulationEffect',
                -rc_name => 'default',
                -max_retry_count => 0,
                -hive_capacity  =>  50,
            }, 
            {   -logic_name => 'finish_regulation_effect',
                -module => 'Bio::EnsEMBL::Variation::Pipeline::FinishRegulationEffect',
                -rc_name => 'long',
                -hive_capacity => 1,
                -rc_name => 'long',
                -parameters => {
                  'pipeline_dir' => $self->o('pipeline_dir'),
                },
            },
        );
    }
    return \@analyses;
}

1;
