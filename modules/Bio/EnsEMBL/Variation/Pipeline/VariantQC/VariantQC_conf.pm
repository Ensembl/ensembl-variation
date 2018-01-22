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



=head1 NAME 

Bio::EnsEMBL::Variation::Pipeline::VariantQC::VariantQC_conf

=head1 DESCRIPTION

Configuration module for variant QC eHive process

=cut

package Bio::EnsEMBL::Variation::Pipeline::VariantQC::VariantQC_conf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');

sub default_options {
    my ($self) = @_;

    # the hash returned from this function is used to configure the pipeline, you can supply
    # any of these options on the command line to override these default values
    
    # you shouldn't need to edit anything in this file other than these values, if you
    # find you do need to then we should probably make it an option here, contact
    # the variation team to discuss this - patches are welcome!

    my $login = `whoami`;
    chomp $login;

    return {

        # general pipeline options that you should change to suit your environment

        hive_use_triggers       => 0,  
        hive_force_init         => 1,
        hive_use_param_stack    => 0,
        hive_auto_rebalance_semaphores => 0,
        hive_no_init            => 0,


        compile_module_once     => 1, 

        ## check why jobs failed before re-submitting for analysis
        retry_throwing_job      => 0,  

        # the location of your checkout of the ensembl API (the hive looks for SQL files here)
        
        ensembl_cvs_root_dir    => $ENV{'HOME'}.'/bin/',
        hive_root_dir           => $ENV{'HOME'}.'/bin/ensembl-hive', 

        # a name for your pipeline (will also be used in the name of the hive database)
        
        pipeline_name           => 'variation_qc',

        # a directory to keep hive output files and your registry file, you should
        # create this if it doesn't exist

        pipeline_dir            => '/gpfs/nobackup/ensembl/' . $login . '/'.$self->o('pipeline_name') . '/'.  $self->o('species'),
#        pipeline_dir => '/hps/nobackup/production/ensembl/anja/release_90/human/37/variant_qc/',
        # a directory where hive workers will dump STDOUT and STDERR for their jobs
        # if you use lots of workers this directory can get quite big, so it's
        # a good idea to keep it on lustre, or some other place where you have a 
        # healthy quota!
        
        output_dir              => $self->o('pipeline_dir').'/hive_output',

        # a standard ensembl registry file containing connection parameters
        # for your target database(s) (and also possibly aliases for your species
        # of interest that you can then supply to init_pipeline.pl with the -species
        # option)
        
        reg_file                => $self->o('pipeline_dir').'/ensembl.registry',


        ## number of *variants* handled per batch

        qc_batch_size            => 1000, 
        unmapped_batch_size      => 100000, ## quicker check can be binned in bigger chunks


        # Options to change for failure recovery 

        ## only data with variation_id >= start_at_variation_id will be imported
                            
        start_at_variation_id     => 1,

        ## this can be changed for failure recovery 
        ## working tables will not be created when create_working_table is set to 0
                                    
        create_working_tables     => 1,

        # create tmp_map_weight table unless this set to 0

        create_map_table          => 1,  


        # configuration for the various resource options used in the pipeline
        # EBI farm users should either change these here, or override them on the
        # command line to suit the EBI farm. The names of each option hopefully
        # reflect their usage, but you may want to change the details (memory
        # requirements, queue parameters etc.) to suit your own data
        
        default_lsf_options => '-R"select[mem>2000] rusage[mem=2000]" -M2000',       ## upped to 4G from 2G for mouse
        urgent_lsf_options  => '-R"select[mem>2000] rusage[mem=2000]" -M2000',
        highmem_lsf_options => '-R"select[mem>15000] rusage[mem=15000]" -M15000', 
        long_lsf_options    => '-R"select[mem>2000] rusage[mem=2000]" -M2000',
        medium_lsf_options  => '-R"select[mem>5000] rusage[mem=5000]" -M5000', ## switched from 4->5 on moving to farm3

        # options controlling the number of workers used for the parallelisable analyses

        variant_qc_capacity        => 60,
        unmapped_var_capacity      => 10,


        # these flags control which parts of the pipeline are run

        run_check_dbSNP_import           => 1, 
        run_create_seqdb                 => 1, 

        run_variant_qc                   => 1,
        run_unmapped_var                 => 1,

        run_flip_population_genotype     => 1,
        run_update_population_genotype   => 1,

        run_PAR_check                    => 1, 
        run_Pubmed_check                 => 1,

        run_evidence_check               => 1,


        # put back support for re-runs on new format schema

        schema                       => 'old',

        # connection parameters for the hive database, you should supply the hive_db_password
        # option on the command line to init_pipeline.pl (parameters for the target database
        # should be set in the registry file defined above)

        # init_pipeline.pl will create the hive database on this machine, naming it
        # <username>_<pipeline_name>, and will drop any existing database with this
        # name

        hive_db_host    => 'mysql-ens-var-prod-2',
        hive_db_port    => 4521,
        hive_db_user    => 'ensadmin',

        pipeline_db => {
            -host   => $self->o('hive_db_host'),
            -port   => $self->o('hive_db_port'),
            -user   => $self->o('hive_db_user'),
            -pass   => $self->o('hive_db_password'),            
            -dbname => $ENV{'USER'}.'_'.$self->o('pipeline_name') . '_' . $self->o('species'),
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
        'medium'  => { 'LSF' => $self->o('medium_lsf_options')  },
    };
}

sub pipeline_analyses {
    my ($self) = @_;

    my @common_params = (
        ensembl_registry    => $self->o('reg_file'),
        species             => $self->o('species'),
        pipeline_dir        => $self->o('pipeline_dir'),
    );
   
   
    my @analyses;

 
    push @analyses, (
    
      { -logic_name => 'init_run_variant_qc',
        -module     => 'Bio::EnsEMBL::Variation::Pipeline::VariantQC::InitVariantQC',
        -parameters => {
            qc_batch_size                  => $self->o('qc_batch_size'),
            unmapped_batch_size            => $self->o('unmapped_batch_size'),

            run_create_seqdb               => $self->o('run_create_seqdb'),
            run_check_dbSNP_import         => $self->o('run_check_dbSNP_import'),
            run_variant_qc                 => $self->o('run_variant_qc'),
            run_unmapped_var               => $self->o('run_unmapped_var'),
            run_flip_population_genotype   => $self->o('run_flip_population_genotype'),
            run_update_population_genotype => $self->o('run_update_population_genotype'),

            start_at_variation_id          => $self->o('start_at_variation_id'),
            create_working_tables          => $self->o('create_working_tables'),
            create_map_table               => $self->o('create_map_table'),
            @common_params,
        },
        -input_ids  => [{}],
        -hive_capacity  => -1,
        -rc_name    => 'default',
        -flow_into  => {
             2 => [ 'check_dbSNP_import' ],
             3 => [ 'create_seqdb' ],
             4 => [ 'variant_qc'     ],
             5 => [ 'unmapped_var'   ],
             6 => [ 'flip_population_genotype' ],
             7 => [ 'update_population_genotype' ],
             8 => [ 'special_cases' ],
             9 => [ 'finish_variation_qc' ],
                },
     },

     {  -logic_name => 'check_dbSNP_import',
        -module     => 'Bio::EnsEMBL::Variation::Pipeline::VariantQC::CheckdbSNPImport',
        -parameters => {    
            @common_params,
        },
        -input_ids      => [],
        -hive_capacity  => -1,
        -rc_name        => 'default',               
      },

     {  -logic_name => 'create_seqdb',
        -module     => 'Bio::EnsEMBL::Variation::Pipeline::VariantQC::CreateSeqDB',
        -parameters => {    
            @common_params,
        },
        -input_ids      => [],
        -hive_capacity  => 1,
        -rc_name        => 'medium',               
      },

          
      {  -logic_name     => 'unmapped_var',
         -module         => 'Bio::EnsEMBL::Variation::Pipeline::VariantQC::UnmappedVariant',        
         -parameters     => {   
               batch_size => $self->o('unmapped_batch_size'),
               @common_params,
           },
         -input_ids        => [],
         -hive_capacity    => $self->o('unmapped_var_capacity'),
#         -analysis_capacity=>$self->o('unmapped_var_capacity'),
         -max_retry_count  => 0,
         -rc_name          => 'default',
         -wait_for         => [ 'check_dbSNP_import' ],
         -flow_into        => {},
         },



     {   -logic_name     => 'variant_qc',
         -module         => 'Bio::EnsEMBL::Variation::Pipeline::VariantQC::VariantQC',
         -parameters     => {   
             schema             => $self->o('schema'),
             batch_size         => $self->o('qc_batch_size'),
             use_seqdb          => $self->o('run_create_seqdb'),  
             evidence_check     => $self->o('run_evidence_check'),  
             @common_params,
         },
         -input_ids        => [],
         -hive_capacity    => $self->o('variant_qc_capacity') ,  ## switch this off on hive upgrade
#         -analysis_capacity=> $self->o('variant_qc_capacity') ,
         -max_retry_count  => 0,
         -rc_name          => 'default',
         -wait_for         => [ 'check_dbSNP_import', 'create_seqdb'],
         -flow_into        => {},
         },
        

    {   -logic_name     => 'flip_population_genotype',
        -module         => 'Bio::EnsEMBL::Variation::Pipeline::VariantQC::FlipPopulationGenotype',
        -parameters     => {                   
            @common_params,
        },
        -input_ids      => [],
        -hive_capacity  => -1,
        -rc_name        => 'long',
        -wait_for       => [  'variant_qc' ],
        -flow_into      => {},
    },



    {   -logic_name     => 'update_population_genotype',
        -module         => 'Bio::EnsEMBL::Variation::Pipeline::VariantQC::UpdatePopulationGenotype', 
        -parameters     => {
            @common_params,
        },
        -input_ids      => [],
        -hive_capacity  => -1,
        -rc_name        => 'default',
        -wait_for       => [ 'variant_qc', 'flip_population_genotype' ],
        -flow_into      => {},
    },              


   {     -logic_name     => 'special_cases',
         -module         => 'Bio::EnsEMBL::Variation::Pipeline::VariantQC::SpecialCase',
         -parameters     => {   
             run_PAR_check      => $self->o('run_PAR_check'),
             run_Pubmed_check   => $self->o('run_Pubmed_check'),
             @common_params,
         },
         -input_ids      => [],
         -hive_capacity  => -1,
         -rc_name        => 'default',
         -wait_for       => [ 'variant_qc'],
         -flow_into      => {},
         },


    {   -logic_name     => 'finish_variation_qc',
        -module         => 'Bio::EnsEMBL::Variation::Pipeline::VariantQC::FinishVariantQC', 
        -parameters     => {
            @common_params,
         },
         -input_ids      => [],
         -hive_capacity  => -1,
         -rc_name        => 'default',
         -wait_for       => [ 'variant_qc','unmapped_var','update_population_genotype','special_cases' ],
         -flow_into      => {},
    },                   
   
    );
    

   
    return \@analyses;
}

1;

