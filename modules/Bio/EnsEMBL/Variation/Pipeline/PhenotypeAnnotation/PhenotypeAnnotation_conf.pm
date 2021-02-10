=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2021] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     https://www.apache.org/licenses/LICENSE-2.0

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

package Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::PhenotypeAnnotation_conf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');

use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
use Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::Constants qw(RGD ANIMALQTL ZFIN GWAS OMIA EGA ORPHANET MIMMORBID DDG2P CGC IMPC MGI NONE HUMAN MOUSE ANIMALSET);


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

        hive_force_init => 1,
        hive_use_param_stack => 0,
        hive_use_triggers => 0,
        hive_auto_rebalance_semaphores => 0,
        hive_no_init => 0,
        hive_debug_init                => 0,     # setting it to 1 will make init_pipeline.pl tell everything it's doing
        hive_default_max_retry_count   => 3,     # default value for the max_retry_count parameter of each analysis

        debug_mode   => 0,

        # the location of your checkout of the ensembl API (the hive looks for SQL files here)

        ensembl_cvs_root_dir    => $ENV{'HOME'} . '/src',
        hive_root_dir           => $ENV{'HOME'} . '/src/ensembl-hive',

        # release number used in the name of default workdir and hive database
        ensembl_release         => 95,

        # a name for your pipeline (will also be used in the name of the hive database)
        
        pipeline_name           => 'phenotype_annotation',

        # a directory to keep hive output files and your registry file, you should
        # create this if it doesn't exist

        pipeline_dir            => '/hps/nobackup2/production/ensembl/' . $login . '/' . $self->o('pipeline_name')."_".$self->o('ensembl_release') ,

        # a directory where hive workers will dump STDOUT and STDERR for their jobs
        # if you use lots of workers this directory can get quite big, so it's
        # a good idea to keep it on lustre, or some other place where you have a
        # healthy quota!

        output_dir              => $self->o('pipeline_dir').'/hive_output',

        # a standard ensembl registry file containing connection parameters
        # for your target database(s) (and also possibly aliases for your species
        # of interest)

        reg_file                => $self->o('pipeline_dir').'/ensembl.registry',

        # the run type can be one of: RGD (import RGD data),
        # AnimalQTL (import AnimalQTL), ZFIN (import ZFIN data)
        # The species which are imported for each data sources are in Constants.pm

        run_type                =>  NONE,

        threshold_qtl           =>  undef, # default for RGD_qtl, AnimalQTL

        ega_database_conf       => $self->o('pipeline_dir').'/ega_database.conf',

        # configuration for the various resource options used in the pipeline
        # Users of other farms should change these here, or override them on
        # the command line to suit your farm. The names of each option hopefully
        # reflect their usage, but you may want to change the details (memory
        # requirements, queue parameters etc.) to suit your own data

        default_lsf_options => '-qproduction-rh74 -R"select[mem>2000] rusage[mem=2000]" -M2000',
        medmem_lsf_options  => '-qproduction-rh74 -R"select[mem>4000] rusage[mem=4000]" -M4000',
        urgent_lsf_options  => '-qproduction-rh74 -R"select[mem>2000] rusage[mem=2000]" -M2000',
        highmem_lsf_options => '-qproduction-rh74 -R"select[mem>15000] rusage[mem=15000] span[hosts=1]" -M15000 -n4', # this is LSF speak for "give me 15GB of memory"
        long_lsf_options    => '-qproduction-rh74 -R"select[mem>2000] rusage[mem=2000]" -M2000',

        # connection parameters for the hive database, you should supply the hive_db_password
        # option on the command line to init_pipeline.pl (parameters for the target database
        # should be set in the registry file defined above)

        hive_use_triggers       => 0,

        # init_pipeline.pl will create the hive database on this machine, naming it
        # <username>_<pipeline_name>_<ensembl_release>, and will drop any existing database with this
        # name

        hive_db_host    => 'mysql-ens-var-prod-1',
        hive_db_port    => 4449,
        hive_db_user    => 'ensadmin',
        hive_db_name    => $ENV{'USER'}.'_ehive_'.$self->o('pipeline_name').'_'.$self->o('ensembl_release'),

        pipeline_db => {
            -host   => $self->o('hive_db_host'),
            -port   => $self->o('hive_db_port'),
            -user   => $self->o('hive_db_user'),
            -pass   => $self->o('hive_db_password'),
            -dbname => $self->o('hive_db_name'),
            -driver => 'mysql',
        },
    };
}


sub resource_classes {
    my ($self) = @_;
    return {
          'default' => { 'LSF' => $self->o('default_lsf_options') },
          'urgent'  => { 'LSF' => $self->o('urgent_lsf_options')  },
          'highmem' => { 'LSF' => $self->o('highmem_lsf_options') },
          'long'    => { 'LSF' => $self->o('long_lsf_options')    },
          'medmem'  => { 'LSF' => $self->o('medmem_lsf_options') },
    };
}

sub pipeline_analyses {
    my ($self) = @_;

    my @common_params = (
        ensembl_registry    => $self->o('reg_file'),
        pipeline_dir        => $self->o('pipeline_dir'),
        debug_mode          => $self->o('debug_mode'),
    );

    return [
        {   -logic_name => 'init_import_phenotype',
            -module     => 'Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::InitPhenotypeAnnotation',
            -parameters => {
                @common_params,
                run_type => $self->o('run_type'),
            },
            -input_ids  => [], #default
            -rc_name    => 'default',
            -max_retry_count => 0,
            -flow_into  => {
                '2' => [ 'import_human' ],
                '3' => [ 'import_mouse' ],
                '4' => [ 'import_animalset' ],

                '5' => [ 'import_rgd' ],
                '6' => [ 'import_zfin' ],
            },
        },


        # HUMAN import:
        {   -logic_name => 'import_human',
            -module     => 'Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::ImportHuman',
            -parameters => {
                @common_params,
            },
            -input_ids      => [], #default
            -hive_capacity  => 1,
            -rc_name    => 'default',
            -flow_into  => {
                '2->A' => [ 'import_gwas' ],
                '3->A' => [ 'import_ega' ],
                '4->A' => [ 'import_orphanet' ],
                '5->A' => [ 'import_mimmorbid' ],
                '6->A' => [ 'import_ddg2p' ],
                '7->A' => [ 'import_cancerGC' ],
                'A->8' => [ 'check_phenotypes'],
            },
        },

        {   -logic_name => 'import_gwas',
            -module     => 'Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::ImportGWAS',
            -parameters => {
                @common_params,
            },
            -input_ids      => [], #default
            -hive_capacity  => 1,
            -rc_name    => 'default',
            -flow_into  => {
                2 => [ 'check_gwas']
            },
        },

        {   -logic_name => 'check_gwas',
            -module     => 'Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::CheckPhenotypeAnnotation',
            -parameters => {
                @common_params,
            },
            -input_ids      => [], #default
            -hive_capacity  => 1,
            -rc_name    => 'default',
            -flow_into  => {
                2 => [ 'import_ega']
            },
            -max_retry_count => 0,
        },

        {   -logic_name => 'import_ega',
            -module     => 'Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::ImportEGA',
            -parameters => {
                ega_database_conf => $self->o('ega_database_conf'),
                @common_params,
            },
            -input_ids      => [], #default
            -hive_capacity  => 1,
            -rc_name    => 'default',
            -flow_into  => {
                2 => [ 'check_ega']
            },
        },

        {   -logic_name => 'check_ega',
            -module     => 'Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::CheckPhenotypeAnnotation',
            -parameters => {
                @common_params,
            },
            -input_ids      => [], #default
            -hive_capacity  => 1,
            -rc_name    => 'default',
            -flow_into  => {
                2 => [ 'import_orphanet']
            },
            -max_retry_count => 0,
        },

        {   -logic_name => 'import_orphanet',
            -module     => 'Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::ImportOrphanet',
            -parameters => {
                @common_params,
            },
            -input_ids      => [], #default
            -hive_capacity  => 1,
            -rc_name    => 'default',
            -flow_into  => {
                2 => [ 'check_orphanet']
            },
        },

        {   -logic_name => 'check_orphanet',
            -module     => 'Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::CheckPhenotypeAnnotation',
            -parameters => {
                @common_params,
            },
            -input_ids      => [], #default
            -hive_capacity  => 1,
            -rc_name    => 'default',
            -flow_into  => {
                2 => [ 'import_mimmorbid']
            },
            -max_retry_count => 0,
        },

        {   -logic_name => 'import_mimmorbid',
            -module     => 'Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::ImportMIMmorbid',
            -parameters => {
                @common_params,
            },
            -input_ids      => [], #default
            -hive_capacity  => 1,
            -rc_name    => 'default',
            -flow_into  => {
                2 => [ 'check_mimmorbid']
            },
        },

        {   -logic_name => 'check_mimmorbid',
            -module     => 'Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::CheckPhenotypeAnnotation',
            -parameters => {
                @common_params,
            },
            -input_ids      => [], #default
            -hive_capacity  => 1,
            -rc_name    => 'default',
            -flow_into  => {
                2 => [ 'import_ddg2p']
            },
            -max_retry_count => 0,
        },

        {   -logic_name => 'import_ddg2p',
            -module     => 'Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::ImportDDG2P',
            -parameters => {
                @common_params,
            },
            -input_ids      => [], #default
            -hive_capacity  => 1,
            -rc_name    => 'default',
            -flow_into  => {
                2 => [ 'check_ddg2p']
            },
        },

        {   -logic_name => 'check_ddg2p',
            -module     => 'Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::CheckPhenotypeAnnotation',
            -parameters => {
                @common_params,
            },
            -input_ids      => [], #default
            -hive_capacity  => 1,
            -rc_name    => 'default',
            -flow_into  => {
                2 => [ 'import_cancerGC']
            },
            -max_retry_count => 0,
        },

        {   -logic_name => 'import_cancerGC',
            -module     => 'Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::ImportCGC',
            -parameters => {
                @common_params,
            },
            -input_ids      => [], #default
            -hive_capacity  => 1,
            -rc_name    => 'default',
            -flow_into  => {
                2 => [ 'check_cancerGC']
            },
        },

        {   -logic_name => 'check_cancerGC',
            -module     => 'Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::CheckPhenotypeAnnotation',
            -parameters => {
                @common_params,
            },
            -input_ids      => [], #default
            -hive_capacity  => 1,
            -rc_name    => 'default',
            -max_retry_count => 0,
        },


        # Mouse import:
        {   -logic_name => 'import_mouse',
            -module     => 'Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::ImportMouse',
            -parameters => {
                @common_params,
            },
            -input_ids      => [], #default
            -hive_capacity  => 1,
            -rc_name    => 'default',
            -flow_into  => {
                '2->A' => [ 'import_impc'],
                '3->A' => [ 'import_mgi'],
                'A->4' => [ 'check_phenotypes'],
            },
        },

        {   -logic_name => 'import_impc',
            -module     => 'Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::ImportIMPC',
            -parameters => {
                @common_params,
            },
            -input_ids      => [], #default
            -hive_capacity  => 1,
            -rc_name    => 'default',
            -flow_into  => {
                2 => { 'check_impc' => INPUT_PLUS() },
            },
        },

        {   -logic_name => 'check_impc',
            -module     => 'Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::CheckPhenotypeAnnotation',
            -parameters => {
                @common_params,
            },
            -input_ids      => [], #default
            -hive_capacity  => 1,
            -rc_name    => 'default',
            -flow_into  => {
                2 =>  { 'import_mgi' => INPUT_PLUS() },
            },
            -max_retry_count => 0,
        },

        {   -logic_name => 'import_mgi',
            -module     => 'Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::ImportMGI',
            -parameters => {
                @common_params,
            },
            -input_ids      => [], #default
            -hive_capacity  => 1,
            -rc_name    => 'default',
            -flow_into  => {
                2 => [ 'check_mgi'],
            },
        },

        {   -logic_name => 'check_mgi',
            -module     => 'Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::CheckPhenotypeAnnotation',
            -parameters => {
                @common_params,
            },
            -input_ids      => [], #default
            -hive_capacity  => 1,
            -rc_name    => 'default',
            -max_retry_count => 0,
        },


        # AnimalSet import:
        {   -logic_name => 'import_animalset',
            -module     => 'Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::ImportAnimalSet',
            -parameters => {
                @common_params,
            },
            -input_ids      => [], #default
            -hive_capacity  => 1,
            -rc_name    => 'default',
            -flow_into  => {
                '2' => [ 'import_omia' ],
                '3' => [ 'import_animalqtldb' ],
            },
        },

        {   -logic_name => 'import_omia',
            -module     => 'Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::ImportOMIA',
            -parameters => {
                @common_params,
            },
            -input_ids      => [], #default
            -hive_capacity  => 1,
            -analysis_capacity => 1,
            -rc_name    => 'default',
            -flow_into  => {
                2 => [ 'check_omia']
            },
        },

        {   -logic_name => 'check_omia',
            -module     => 'Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::CheckPhenotypeAnnotation',
            -parameters => {
                @common_params,
            },
            -input_ids      => [], #default
            -hive_capacity  => 1,
            -analysis_capacity => 1,
            -rc_name    => 'default',
            -flow_into  => {
                2 => [ 'import_animalqtldb'],
                3 => [ 'check_phenotypes'],
            },
            -max_retry_count => 0,
        },

        {   -logic_name => 'import_animalqtldb',
            -module     => 'Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::ImportAnimalQTL',
            -parameters => {
                threshold_qtl       => $self->o('threshold_qtl'),
                @common_params,
            },
            -input_ids      => [], #default
            -hive_capacity  => 1,
            -analysis_capacity => 1,
            -rc_name    => 'default',
            -flow_into  => {
                2 => [ 'check_animalqtl']
            },
        },

        {   -logic_name => 'check_animalqtl',
            -module     => 'Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::CheckPhenotypeAnnotation',
            -parameters => {
                @common_params,
            },
            -input_ids      => [], #default
            -hive_capacity  => 1,
            -analysis_capacity => 1,
            -rc_name    => 'default',
            -flow_into => {
                2 => [ 'check_phenotypes']
            },
            -max_retry_count => 0,
        },


        #Other species:
        {   -logic_name => 'import_rgd',
            -module     => 'Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::ImportRGD',
            -parameters => {
                threshold_qtl   => $self->o('threshold_qtl'),
                @common_params,
            },
            -input_ids      => [], #default
            -hive_capacity  => 1,
            -rc_name    => 'default',
            -flow_into  => {
                2 => [ 'check_phenotypes']
            },
        },

        {   -logic_name => 'import_zfin',
            -module     => 'Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::ImportZFIN',
            -parameters => {
                @common_params,
            },
            -input_ids      => [], #default
            -hive_capacity  => 1,
            -rc_name    => 'default',
            -flow_into  => {
                2 => [ 'check_phenotypes']
            },
        },

        {   -logic_name => 'check_phenotypes',
            -module     => 'Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::CheckPhenotypeAnnotation',
            -parameters => {
                @common_params,
            },
            -input_ids      => [], #default
            -hive_capacity  => 1,
            -rc_name    => 'default',
            -flow_into  => {
                2 => [ 'import_ontology_mapping'],
                3 => [ 'finish_phenotype_annotation']
            },
        },

        {   -logic_name => 'import_ontology_mapping',
            -module     => 'Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::OntologyMapping',
            -parameters => {
                @common_params,
            },
            -input_ids      => [], #default
            -hive_capacity  => 1,
            -analysis_capacity => 1,
            -rc_name    => 'default',
            -flow_into  => {
                2 => [ 'finish_phenotype_annotation']
            },
        },

        {   -logic_name => 'finish_phenotype_annotation',
            -module     => 'Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::FinishPhenotypeAnnotation',
            -parameters => {
                @common_params,
            },
            -input_ids      => [], #default
            -hive_capacity  => 1,
            -analysis_capacity => 1,
            -rc_name    => 'default',
            -flow_into      => {},
            -failed_job_tolerance => 0,
            -max_retry_count => 0,
        },
        
        {   -logic_name => 'finish_pipeline',
            -module     => 'Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::FinishPhenotypeAnnotationPipeline',
            -parameters => {
                pipeline_name => $self->o('pipeline_name'),
                @common_params,
            },
            -input_ids      => [], #default
            -hive_capacity  => 1,
            -rc_name    => 'default',
            -flow_into      => {}, #default
            -failed_job_tolerance => 0,
            -max_retry_count => 0,
        },
      ];
}

1;
