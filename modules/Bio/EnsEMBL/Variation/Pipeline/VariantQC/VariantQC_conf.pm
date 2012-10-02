=head1 LICENSE

 Copyright (c) 1999-2012 The European Bioinformatics Institute and
 Genome Research Limited.  All rights reserved.

 This software is distributed under a modified Apache license.
 For license details, please see

   http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <dev@ensembl.org>.

 Questions may also be sent to the Ensembl help desk at
 <helpdesk@ensembl.org>.



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

        ## check why jobs failed before re-submitting for analysis
        retry_throwing_job      => 0,  

        # the location of your checkout of the ensembl API (the hive looks for SQL files here)
        
        ensembl_cvs_root_dir    => $ENV{'HOME'}.'/EBI/bin/HEAD',

        # a name for your pipeline (will also be used in the name of the hive database)
        
        pipeline_name           => 'variation_qc',

        # a directory to keep hive output files and your registry file, you should
        # create this if it doesn't exist

        pipeline_dir            => '/lustre/scratch110/ensembl/' . $login . '/'.$self->o('pipeline_name') . '/'.  $self->o('species'),

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

        # flip varation_features on reverse strand with map weight 1 unless this set to 0

        #flip_variants            => 1,    TO BE IMPLEMENTED on by default

        # The current dbSNP importer does not create the variation_feature.allele_string
        # Switch this to 0 if QC'ing external data imported with variation_feature.allele_string's

        do_allele_string          => 1,

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
        
        default_lsf_options => '-R"select[mem>2000] rusage[mem=2000]" -M2000000',
        urgent_lsf_options  => '-R"select[mem>2000] rusage[mem=2000]" -M2000000',
        highmem_lsf_options => '-R"select[mem>15000] rusage[mem=15000]" -M15000000', # this is Sanger LSF speak for "give me 15GB of memory"
        long_lsf_options    => '-q long -R"select[mem>2000] rusage[mem=2000]" -M2000000',

        # options controlling the number of workers used for the parallelisable analyses

        variant_qc_capacity        => 20,
        unmapped_var_capacity      => 10,


        # these flags control which parts of the pipeline are run

        run_check_dbSNP_import           => 1,
        run_variant_qc                   => 1,
        run_unmapped_var                 => 1,
        run_flip_population_genotype     => 1,
        run_update_population_genotype   => 1,

        # put back support for re-runs on new format schema

        schema                       => 'old',

        # connection parameters for the hive database, you should supply the hive_db_password
        # option on the command line to init_pipeline.pl (parameters for the target database
        # should be set in the registry file defined above)

        # init_pipeline.pl will create the hive database on this machine, naming it
        # <username>_<pipeline_name>, and will drop any existing database with this
        # name

        hive_db_host    => 'ens-variation',
        hive_db_port    => 3306,
        hive_db_user    => 'ensadmin',

        pipeline_db => {
            -host   => $self->o('hive_db_host'),
            -port   => $self->o('hive_db_port'),
            -user   => $self->o('hive_db_user'),
            -pass   => $self->o('hive_db_password'),            
            -dbname => $ENV{'USER'}.'_'.$self->o('pipeline_name') . '_' . $self->o('species'),
        },
    };
}

sub pipeline_create_commands {
    my ($self) = @_;
    return [
        'mysql '.$self->dbconn_2_mysql('pipeline_db', 0).q{-e 'DROP DATABASE IF EXISTS }.$self->o('pipeline_db', '-dbname').q{'},
        @{$self->SUPER::pipeline_create_commands}, 
        'mysql '.$self->dbconn_2_mysql('pipeline_db', 1).q{-e 'INSERT INTO meta (meta_key, meta_value) VALUES ("hive_output_dir", "}.$self->o('output_dir').q{")'},
    ];
}

sub resource_classes {
    my ($self) = @_;
    return {
        'default' => { 'LSF' => $self->o('default_lsf_options') },
        'urgent'  => { 'LSF' => $self->o('urgent_lsf_options')  },
        'highmem' => { 'LSF' => $self->o('highmem_lsf_options') },
        'long'    => { 'LSF' => $self->o('long_lsf_options')    },
    };
}

sub pipeline_analyses {
    my ($self) = @_;

    my @common_params = (
        ensembl_registry    => $self->o('reg_file'),
        species             => $self->o('species'),
    );
   
    my @analyses;

    

    push @analyses, (
    
      { -logic_name => 'init_run_variant_qc',
        -module     => 'Bio::EnsEMBL::Variation::Pipeline::VariantQC::InitVariantQC',
        -parameters => {
            qc_batch_size                  => $self->o('qc_batch_size'),
            unmapped_batch_size            => $self->o('unmapped_batch_size'),

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
        -rc_id      => 'default',
        -flow_into  => {
             2 => [ 'check_dbSNP_import' ],
             3 => [ 'variant_qc'     ],
             4 => [ 'unmapped_var'   ],
             5 => [ 'flip_population_genotype' ],
             6 => [ 'update_population_genotype' ],
             7 => [ 'finish_variation_qc' ],
                },
     },

     {  -logic_name => 'check_dbSNP_import',
        -module     => 'Bio::EnsEMBL::Variation::Pipeline::VariantQC::CheckdbSNPImport',
        -parameters => {    
            pipeline_dir  => $self->o('pipeline_dir'),       
            @common_params,
        },
        -input_ids      => [],
        -hive_capacity  => 1,
        -rc_id          => 'default',               
      },

          
      {  -logic_name     => 'unmapped_var',
         -module         => 'Bio::EnsEMBL::Variation::Pipeline::VariantQC::UnmappedVariant',        
         -parameters     => {   
               batch_size => $self->o('unmapped_batch_size'),
               @common_params,
           },
         -input_ids      => [],
         -hive_capacity  => $self->o('unmapped_var_capacity'),
         -rc_id          => 'default',
         -wait_for       => [ 'check_dbSNP_import' ],
         -flow_into      => {},
         },



     {   -logic_name     => 'variant_qc',
         -module         => 'Bio::EnsEMBL::Variation::Pipeline::VariantQC::VariantQC',
         -parameters     => {   
             schema             => $self->o('schema'),
             batch_size         => $self->o('qc_batch_size'),
             do_allele_string   => $self->o('do_allele_string'),
             pipeline_dir       => $self->o('pipeline_dir'),   ### temp - write temp log files
             @common_params,
         },
         -input_ids      => [],
         -hive_capacity  => $self->o('variant_qc_capacity'),
         -rc_id          => 'default',
         -wait_for       => [ 'check_dbSNP_import' ],
         -flow_into      => {},
     },



    {   -logic_name     => 'flip_population_genotype',
        -module         => 'Bio::EnsEMBL::Variation::Pipeline::VariantQC::FlipPopulationGenotype',
        -parameters     => {                   
            @common_params,
        },
        -input_ids      => [],
        -hive_capacity  => 1,
        -rc_id          => 'default',
        -wait_for       => [ 'check_dbSNP_import' ],
        -flow_into      => {},
    },



    {   -logic_name     => 'update_population_genotype',
        -module         => 'Bio::EnsEMBL::Variation::Pipeline::VariantQC::UpdatePopulationGenotype', 
        -parameters     => {
            @common_params,
        },
        -input_ids      => [],
        -hive_capacity  => 1,
        -rc_id          => 'default',
        -wait_for       => [ 'variant_qc', 'flip_population_genotype' ],
        -flow_into      => {},
    },              



    {   -logic_name     => 'finish_variation_qc',
        -module         => 'Bio::EnsEMBL::Variation::Pipeline::VariantQC::FinishVariantQC', 
        -parameters     => {
            pipeline_dir  => $self->o('pipeline_dir'),
            @common_params,
         },
         -input_ids      => [],
         -hive_capacity  => 1,
         -rc_id          => 'default',
         -wait_for       => [ 'variant_qc','unmapped_var','update_population_genotype' ],
         -flow_into      => {},
    },                   
   
    );
    

   
    return \@analyses;
}

1;

