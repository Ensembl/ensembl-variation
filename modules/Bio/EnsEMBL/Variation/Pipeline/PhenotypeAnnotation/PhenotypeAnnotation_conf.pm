=head1 LICENSE

Copyright [2018] EMBL-European Bioinformatics Institute

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
use Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::Constants qw(RGD AnimalQTL NONE);

#TODO: Q: should I use EnsemblGeneric_conf.pm? / clear understanding when to use EnsemblGeneric_conf, when HiveGeneric_conf and when neither/smth else?


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

        debug_mode   => 0,

        # the location of your checkout of the ensembl API (the hive looks for SQL files here)

        ensembl_cvs_root_dir    => $ENV{'HOME'} . '/src',
        hive_root_dir           => $ENV{'HOME'} . '/src/ensembl-hive',

        # release number used in the name of default workdir and hive database
        ensembl_release         => 93,

        # a name for your pipeline (will also be used in the name of the hive database)
        
        pipeline_name           => 'phenotype_annotation',

        # a directory to keep hive output files and your registry file, you should
        # create this if it doesn't exist

        pipeline_dir            => '/hps/nobackup/production/ensembl/' . $login . '/' . $self->o('pipeline_name')."_".$self->o('ensembl_release') ,

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

        run_import_type         =>  NONE,

        threshold_qtl           =>  0, #for RGD_qtl, AnimalQTL

        animalqtl_input_dir     => $self->o('pipeline_dir').'/AnimalQTL/inputFiles',
        animalqtl_version       => '20180822', #release 36 is 20180822 TODO: confirm there is no computational way to get it

        zfin_version            => '20001020', #13 Sep 2018 #TODO: confirm there is no computational way to get it

        nhgri_version           => '20001020', #TODO: confirm there is no computational way to get it
        # configuration for the various resource options used in the pipeline
        # Users of other farms should change these here, or override them on
        # the command line to suit your farm. The names of each option hopefully
        # reflect their usage, but you may want to change the details (memory
        # requirements, queue parameters etc.) to suit your own data

        default_lsf_options => '-qproduction-rh7 -R"select[mem>2000] rusage[mem=2000]" -M2000',
        medmem_lsf_options  => '-qproduction-rh7 -R"select[mem>4000] rusage[mem=4000]" -M4000',
        urgent_lsf_options  => '-qproduction-rh7 -R"select[mem>2000] rusage[mem=2000]" -M2000',
        highmem_lsf_options => '-qproduction-rh7 -R"select[mem>15000] rusage[mem=15000] span[hosts=1]" -M15000 -n4', # this is LSF speak for "give me 15GB of memory"
        long_lsf_options    => '-qproduction-rh7 -R"select[mem>2000] rusage[mem=2000]" -M2000',

        # connection parameters for the hive database, you should supply the hive_db_password
        # option on the command line to init_pipeline.pl (parameters for the target database
        # should be set in the registry file defined above)

        # Should hive use triggeres?
        hive_use_triggers       => 0,

        # init_pipeline.pl will create the hive database on this machine, naming it
        # <username>_<pipeline_name>, and will drop any existing database with this
        # name

        hive_db_host    => 'mysql-ens-var-prod-1',
        hive_db_port    => 4449,
        hive_db_user    => 'ensadmin',

        pipeline_db => {
            -host   => $self->o('hive_db_host'),
            -port   => $self->o('hive_db_port'),
            -user   => $self->o('hive_db_user'),
            -pass   => $self->o('hive_db_password'),
            -dbname => $ENV{'USER'}.'_'.$self->o('pipeline_name').'_'.$self->o('ensembl_release'),
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
                run_import_type   => $self->o('run_import_type'),
                @common_params,
            },
            -input_ids  => [{}],
            -rc_name    => 'long',
            -max_retry_count => 0,
            -flow_into  => {
                '2->A' => [ 'import_rgd' ],
                '3->A' => [ 'import_animal_qtldb' ],
                '4->A' => [ 'import_zfin' ],
                '5->A' => [ 'import_gwas' ],
                'A->1' => [ 'finish_pipeline' ],
          #      4 => [ 'import_mim_morbid' ],
            #    5 => [ 'import_orphanet' ],
            #    6 => [ 'import_ddg2p' ],
            #    7 => [ 'import_cancer' ],
              #import_gwas, import_omia
              #import_uniprot, import_omim, import_ega, import_giant, import_magic
              #import_zfin, import_goa, import_mgp, import_3i, import_mgi
            },
        },

        {   -logic_name => 'import_rgd',
            -module     => 'Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::ImportRGD',
            -parameters => {
                threshold_qtl   => $self->o('threshold_qtl'),
                @common_params,
            },
            -input_ids      => [],
            -hive_capacity  => 1,
            -rc_name    => 'default',
            -flow_into  => {
                1 => [ 'check_phenotypes']
            },
            -failed_job_tolerance => 1, # tries 1 times to run a job
        },

        {   -logic_name => 'import_animal_qtldb',
            -module     => 'Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::ImportAnimalQTL',
            -parameters => {
                animalqtl_input_dir => $self->o('animalqtl_input_dir'),
                animalqtl_version   => $self->o('animalqtl_version'),
                threshold_qtl       => $self->o('threshold_qtl'),
                @common_params,
            },
            -input_ids      => [],
            -hive_capacity  => 1,
            -rc_name    => 'default',
            -flow_into  => {
                1 => [ 'check_phenotypes']
            },
            -failed_job_tolerance => 5, # tries 5 times to run a job
        },

        {   -logic_name => 'import_zfin',
            -module     => 'Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::ImportZFIN',
            -parameters => {
                zfin_version => $self->o('zfin_version'),
                @common_params,
            },
            -input_ids      => [],
            -hive_capacity  => 1,
            -rc_name    => 'default',
            -flow_into  => {
                1 => [ 'check_phenotypes']
            },
            -failed_job_tolerance => 5, # tries 5 times to run a job
        },

        {   -logic_name => 'import_gwas',
            -module     => 'Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::ImportGWAS',
            -parameters => {
                nhgri_version => $self->o('nhgri_version'),
                @common_params,
            },
            -input_ids      => [],
            -hive_capacity  => 1,
            -rc_name    => 'default',
            -flow_into  => {
                1 => [ 'check_phenotypes']
            },
            -failed_job_tolerance => 5, # tries 5 times to run a job
        },

        {   -logic_name => 'check_phenotypes',
            -module     => 'Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::CheckPhenotypeAnnotation',
            -parameters => {
                @common_params,
            },
            -input_ids      => [],
            -hive_capacity  => 1,
            -rc_name    => 'default',
            -flow_into  => {
                1 => [ 'finish_phenotype_annotation']
          #      2 => [ 'ontology_mapping']
            },
        },
        #next TODO: ontology_mapping

        {   -logic_name => 'finish_phenotype_annotation',
            -module     => 'Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::FinishPhenotypeAnnotation',
            -parameters => {
                @common_params,
            },
            -input_ids      => [],
            -hive_capacity  => 1,
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
            -input_ids      => [],
            -hive_capacity  => 1,
            -rc_name    => 'default',
        #    -wait_for       => [ 'finish_phenotype_annotation' ],
            -flow_into      => {},
            -failed_job_tolerance => 0,
            -max_retry_count => 0,
        },
      ];
}

1;
