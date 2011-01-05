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

        'output_dir'    => '/lustre/scratch103/ensembl/gr5/variation_consequence/hive_output',

        'pipeline_db' => {
            -host   => 'ens-variation',
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
        2 => { -desc => 'highmem',  'LSF' => '-R"select[mem>5000] rusage[mem=5000]" -M5000000'}
    };
}

sub pipeline_analyses {
    my ($self) = @_;
    return [
        {   -logic_name => 'init_jobs',
            -module     => 'Bio::EnsEMBL::Variation::Pipeline::InitJobs',
            -parameters => {},
            -input_ids  => [{
                pph_dir => '/lustre/scratch103/ensembl/gr5/polyphen',
                ensembl_registry => '/lustre/scratch103/ensembl/gr5/variation_consequence/ensembl.registry',
                species => 'Human',
            }],
            -rc_id      => 0,
            -flow_into  => {
                2 => [ 'transcript_effect' ],
            },
        },

        {   -logic_name     => 'transcript_effect',
            -module         => 'Bio::EnsEMBL::Variation::Pipeline::TranscriptEffect',
            -parameters     => {},
            -input_ids      => [],
            -hive_capacity  => 200,
            -rc_id          => 0,
            -flow_into      => {
                #2 => [ 'run_polyphen' ],
            },
        },

#        {   -logic_name     => 'run_polyphen',
#            -module         => 'Bio::EnsEMBL::Variation::Pipeline::RunPolyPhen',
#            -parameters     => {},
#            -input_ids      => [],
#            -hive_capacity  => 100,
#            -rc_id          => 2,
#            -flow_into      => {
#                3   => [ 'run_weka' ],
#            },
#        },
#        
#        {   -logic_name => 'run_weka',
#            -module     => 'Bio::EnsEMBL::Variation::Pipeline::RunWeka',
#            -parameters => {},
#            -input_ids  => [],
#            -rc_id      => 0,
#            -flow_into  => {},
#        },
    ];
}

1;

