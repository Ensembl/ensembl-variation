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
    return {
        'ensembl_cvs_root_dir' => $ENV{'HOME'}.'/workspace',

        'pipeline_name' => 'variation_consequence',

        'output_dir'    => '/lustre/scratch101/ensembl/gr5/variation_consequence/hive_output',
        
        'reg_file'      => '/lustre/scratch101/ensembl/gr5/variation_consequence/ensembl.registry',

        'pipeline_db' => {
            -host   => 'ens-genomics2',
            -port   => 3306,
            -user   => 'ensadmin',
            -pass   => $self->o('password'),            
            -dbname => 'grsr_'.$self->o('pipeline_name'),
        },
    };
}

sub pipeline_create_commands {
    my ($self) = @_;
    return [
        'mysql '.$self->dbconn_2_mysql('pipeline_db', 1).q{-e 'DROP DATABASE IF EXISTS }.$self->o('pipeline_db', '-dbname').q{'},
        @{$self->SUPER::pipeline_create_commands}, 
        'mysql '.$self->dbconn_2_mysql('pipeline_db', 1).q{-e 'INSERT INTO meta (meta_key, meta_value) VALUES ("hive_output_dir", "}.$self->o('output_dir').q{")'},
    ];
}

sub resource_classes {
    my ($self) = @_;
    return {
        0 => { -desc => 'default',  'LSF' => '' },
        1 => { -desc => 'urgent',   'LSF' => '-q yesterday' },
        2 => { -desc => 'highmem',  'LSF' => '-R"select[mem>15000] rusage[mem=15000]" -M15000000'},
        3 => { -desc => 'long',     'LSF' => '-q long' },
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
                1 => [ 'rebuild_consequence_indexes' ],
                2 => [ 'update_variation_feature' ],
                3 => [ 'init_variation_class' ],
                4 => [ 'transcript_effect' ],
            },
        },

        {   -logic_name     => 'transcript_effect',
            -module         => 'Bio::EnsEMBL::Variation::Pipeline::TranscriptEffect',
            -parameters     => {},
            -input_ids      => [],
            -hive_capacity  => 50,
            -rc_id          => 0,
            -flow_into      => {},
        },

        {   -logic_name     => 'rebuild_transcript_variation_indexes',
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
            -wait_for       => [ 'rebuild_transcript_variation_indexes' ],
            -flow_into      => {},
        },
        
        {   -logic_name     => 'init_variation_class',
            -module         => 'Bio::EnsEMBL::Variation::Pipeline::InitVarClass',
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
            -hive_capacity  => 10,
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

