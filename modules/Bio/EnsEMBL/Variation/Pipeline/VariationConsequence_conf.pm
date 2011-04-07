=head1 LICENSE

 Copyright (c) 1999-2011 The European Bioinformatics Institute and
 Genome Research Limited.  All rights reserved.

 This software is distributed under a modified Apache license.
 For license details, please see

   http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <dev@ensembl.org>.

 Questions may also be sent to the Ensembl help desk at
 <helpdesk@ensembl.org>.

=cut

package Bio::EnsEMBL::Variation::Pipeline::VariationConsequence_conf;

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

    return {

        # general pipeline options that you should change to suit your environment
        
        # the location of your checkout of the ensembl API
        
        ensembl_cvs_root_dir    => $ENV{'HOME'}.'/workspace',

        # a name for your pipeline (will also be used in the name of the hive database)
        
        pipeline_name           => 'variation_consequence',

        # a directory where hive workers will dump STDOUT and STDERR for their jobs
        # if you use lots of workers this direcotry can get quite big, so it's
        # a good idea to keep it on lustre, or some other place where you have a 
        # healthy quota!
        
        output_dir              => '/lustre/scratch101/ensembl/gr5/variation_consequence/hive_output',

        # a standard ensembl registry file containing connection parameters
        # for your target database(s) (and also possibly aliases for your species
        # of interest that you can then supply to init_pipeline.pl with the -species
        # option)
        
        reg_file                => '/lustre/scratch101/ensembl/gr5/variation_consequence/ensembl.registry',

        # if set to 1 this option tells the transcript_effect analysis to disambiguate
        # ambiguity codes in single nucleotide alleles, so e.g. an allele string like
        # 'T/M' will be treated as if it were 'T/A/C' (this was a request from ensembl
        # genomes)
        
        disambiguate_single_nucleotide_alleles => 0,

        # configuration for the various resource options used in the pipeline
        # EBI users should either change these here, or override them on the
        # command line to suit the EBI farm. The names if each option hopefully
        # reflect their usage, but you may want to change the details (memory
        # requirements, queue parameters etc.) to suit your own data
        
        default_lsf_options => '',
        urgent_lsf_options  => '-q yesterday',
        highmem_lsf_options => '-R"select[mem>15000] rusage[mem=15000]" -M15000000', # this is Sanger LSF speak for "give me 15GB of memory"
        long_lsf_options    => '-q long',

        # options controlling the number of workers used for the parallelisable analyses
        # these default values work for most species, for human a value of 300 for
        # transcript_effect_capacity seems to work well (though can use a lot of database
        # connections)

        transcript_effect_capacity      => 50,
        set_variation_class_capacity    => 10,

        # connection parameters for the hive database, you should supply the hive_db_pass
        # option on the command line to init_pipeline.pl (parameters for the target database
        # should be set in the registry file defined above)

        # init_pipeline.pl will create the hive database on this machine, naming it
        # <username>_<pipeline_name>, and will drop any existing database with this
        # name

        hive_db_host    => 'ens-genomics2',
        hive_db_port    => 3306,
        hive_db_user    => 'ensadmin',

        pipeline_db => {
            -host   => $self->o('hive_db_host'),
            -port   => $self->o('hive_db_port'),
            -user   => $self->o('hive_db_user'),
            -pass   => $self->o('hive_db_password'),            
            -dbname => $ENV{'USER'}.'_'.$self->o('pipeline_name'),
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
        0 => { -desc => 'default',  'LSF' => $self->o('default_lsf_options') },
        1 => { -desc => 'urgent',   'LSF' => $self->o('urgent_lsf_options')  },
        2 => { -desc => 'highmem',  'LSF' => $self->o('highmem_lsf_options') },
        3 => { -desc => 'long',     'LSF' => $self->o('long_lsf_options')    },
    };
}

sub pipeline_analyses {
    my ($self) = @_;

    return [
        {   -logic_name => 'init_transcript_effect',
            -module     => 'Bio::EnsEMBL::Variation::Pipeline::InitTranscriptEffect',
            -parameters => {},
            -input_ids  => [{
                    ensembl_registry    => $self->o('reg_file'),
                    species             => $self->o('species'),
            }],
            -rc_id      => 1,
            -flow_into  => {
                1 => [ 'rebuild_tv_indexes' ],
                2 => [ 'update_variation_feature' ],
                3 => [ 'init_variation_class' ],
                4 => [ 'transcript_effect' ],
            },
        },

        {   -logic_name     => 'transcript_effect',
            -module         => 'Bio::EnsEMBL::Variation::Pipeline::TranscriptEffect',
            -parameters     => { disambiguate_single_nucleotide_alleles => $self->o('disambiguate_single_nucleotide_alleles') },
            -input_ids      => [],
            -hive_capacity  => $self->o('transcript_effect_capacity'),
            -rc_id          => 0,
            -flow_into      => {},
        },

        {   -logic_name     => 'rebuild_tv_indexes',
            -module         => 'Bio::EnsEMBL::Variation::Pipeline::RebuildIndexes',
            -parameters     => {},
            -input_ids      => [],
            -hive_capacity  => 1,
            -rc_id          => 1,
            -wait_for       => [ 'transcript_effect' ],
            -flow_into      => {},
        },
        
        {   -logic_name     => 'update_variation_feature',
            -module         => 'Bio::EnsEMBL::Variation::Pipeline::UpdateVariationFeature',
            -parameters     => {},
            -input_ids      => [],
            -hive_capacity  => 1,
            -rc_id          => 1,
            -wait_for       => [ 'rebuild_tv_indexes' ],
            -flow_into      => {},
        },
        
        {   -logic_name     => 'init_variation_class',
            -module         => 'Bio::EnsEMBL::Variation::Pipeline::InitVariationClass',
            -parameters     => {num_chunks => 50},
            -input_ids      => [],
            -hive_capacity  => 1,
            -rc_id          => 2,
            -wait_for       => [ 'update_variation_feature' ],
            -flow_into      => {
                1 => [ 'finish_variation_class' ],
                2 => [ 'set_variation_class' ],
            },
        },
        
        {   -logic_name     => 'set_variation_class',
            -module         => 'Bio::EnsEMBL::Variation::Pipeline::SetVariationClass',
            -parameters     => {},
            -input_ids      => [],
            -hive_capacity  => $self->o('set_variation_class_capacity'),
            -rc_id          => 0,
            -flow_into      => {},
        },

        {   -logic_name     => 'finish_variation_class',
            -module         => 'Bio::EnsEMBL::Variation::Pipeline::FinishVariationClass',
            -parameters     => {},
            -input_ids      => [],
            -hive_capacity  => 1,
            -rc_id          => 1,
            -wait_for       => [ 'set_variation_class' ],
            -flow_into      => {},
        },
    ];
}

1;

