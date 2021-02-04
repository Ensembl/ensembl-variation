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

package Bio::EnsEMBL::Variation::Pipeline::AncestralAlleles::AncestralAlleles_conf;

use strict;
use warnings;
use Bio::EnsEMBL::Hive::Version 2.5;
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
    
    my $login = `whoami`;
    chomp $login;

    return {

        # general pipeline options that you should change to suit your environment
        hive_debug_init => 1, # If set to 1, will show the objects (analyses, data-flow rules, etc) that are parsed from the PipeConfig file.
        hive_default_max_retry_count => 0,
        hive_force_init => 1,
        hive_use_param_stack => 0,
        hive_use_triggers => 0,
        hive_auto_rebalance_semaphores => 0, 
        hive_no_init => 0,
        # the location of your checkout of the ensembl API (the hive looks for SQL files here)
        
        ensembl_cvs_root_dir    => $ENV{'HOME'} . '/bin',
        hive_root_dir           => $ENV{'HOME'} . '/bin/ensembl-hive', 
        # a name for your pipeline (will also be used in the name of the hive database)
        
        pipeline_name           => 'ancestral_alleles',

        # a directory to keep hive output files and your registry file, you should
        # create this if it doesn't exist

        pipeline_dir            => $self->o('pipeline_dir'),
        compara_dir             => $self->o('compara_dir'),

        reg_file                => $self->o('pipeline_dir') . '/ensembl.registry',
        batch_size => 1_000_000,
        create_stats => 0,
        non_dbSNP_only => 0,
        
        
        default_lsf_options => '-qproduction-rh74 -R"select[mem>2000] rusage[mem=2000]" -M2000',
        medmem_lsf_options  => '-qproduction-rh74 -R"select[mem>5000] rusage[mem=5000]" -M5000',

        hive_db_host    => 'mysql-ens-var-prod-1',
        hive_db_port    => 4449,
        hive_db_user    => 'ensadmin',

        pipeline_db => {
            -host   => $self->o('hive_db_host'),
            -port   => $self->o('hive_db_port'),
            -user   => $self->o('hive_db_user'),
            -pass   => $self->o('hive_db_password'),            
            -dbname => $ENV{'USER'} . '_ehive_' . $self->o('pipeline_name') . '_' . $self->o('ensembl_release'),
            -driver => 'mysql',
            -reconnect_when_lost => 1
        },
    };
}

sub resource_classes {
    my ($self) = @_;
    return {
          'default' => { 'LSF' => $self->o('default_lsf_options') },
          'medmem'  => { 'LSF' => $self->o('medmem_lsf_options') },
    };
}

sub pipeline_analyses {
    my ($self) = @_;

    my @common_params = (
        ensembl_registry => $self->o('reg_file'),
        pipeline_dir     => $self->o('pipeline_dir'),
        compara_dir      => $self->o('compara_dir'),
        batch_size       => $self->o('batch_size'),
        create_stats     => $self->o('create_stats'),
        non_dbSNP_only   => $self->o('non_dbSNP_only'),

    );
   
    my @analyses;
    push @analyses, (
        {   -logic_name     => 'init',
            -module         => 'Bio::EnsEMBL::Variation::Pipeline::AncestralAlleles::Init',
            -parameters     => {
                @common_params,
            },
            -input_ids      => [{}],
            -rc_name        => 'default',
            -analysis_capacity => 45,
            -hive_capacity => 45,
            -flow_into => { 
              '2->A' => ['assign'],
              'A->1' => ['post_processing']
            },    
        },
        {   -logic_name     => 'assign',
            -module         => 'Bio::EnsEMBL::Variation::Pipeline::AncestralAlleles::Assign',
            -analysis_capacity => 45,
            -hive_capacity => 45,
            -parameters     => {
                @common_params,
            },
            -rc_name        => 'default',
        },
        {   -logic_name     => 'post_processing',
            -module         => 'Bio::EnsEMBL::Variation::Pipeline::AncestralAlleles::PostProcessing',
            -parameters     => {
                @common_params,
            },
            -hive_capacity  => 1,
            -rc_name        => 'default',
        },
    );
    return \@analyses;
}

1;

