package Bio::EnsEMBL::Variation::Pipeline::ProteinFunction::ProteinFunction_conf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');

sub default_options {
    my ($self) = @_;

    return {
    
        ensembl_cvs_root_dir    => $ENV{'HOME'}.'/workspace',

        pipeline_name           => 'protein_function',

        pipeline_dir            => '/lustre/scratch101/ensembl/gr5/protein_function',
        
        output_dir              => $self->o('pipeline_dir').'/hive_output',

        ensembl_registry        => $self->o('pipeline_dir').'/ensembl.registry',

        proteins_fasta          => $self->o('pipeline_dir').'/ensembl_63_proteins.fa',

        species                 => 'human',

        pph_dir                 => $self->o('pipeline_dir').'/polyphen-2.1.0',
        
        sift_dir                => $self->o('pipeline_dir').'/sift-4.0.3',
        
        ncbi_dir                => '/software/ncbiblast/bin',
        
        blastdb                 => $self->o('pipeline_dir').'/blastdb/swissprot_trembl.uni',

        pipeline_db => {
            -host   => 'ens-variation',
            -port   => 3306,
            -user   => 'ensadmin',
            -pass   => $self->o('password'),            
            -dbname => 'grsr_'.$self->o('pipeline_name').'_hive',
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
        0 => { -desc => 'default',  'LSF' => '' },
        1 => { -desc => 'urgent',   'LSF' => '-q yesterday' },
        2 => { -desc => 'highmem',  'LSF' => '-R"select[mem>5000] rusage[mem=5000]" -M5000000 -q long'},
        3 => { -desc => 'long',     'LSF' => '-q long' },
    };
}

sub pipeline_analyses {
    my ($self) = @_;
    
    my @common_params = (
        proteins_fasta      => $self->o('proteins_fasta'),
        ensembl_registry    => $self->o('ensembl_registry'),
        species             => $self->o('species'),
    );

    return [
        {   -logic_name => 'init_jobs',
            -module     => 'Bio::EnsEMBL::Variation::Pipeline::ProteinFunction::InitJobs',
            -parameters => {
                @common_params,
            },
            -input_ids  => [{}],
            -rc_id      => 2,
            -flow_into  => {
                1 => [ 'rebuild_polyphen_indexes' ],
                2 => [ 'run_polyphen' ],
                3 => [ 'rebuild_sift_indexes' ],
                4 => [ 'run_sift' ],
            },
        },

        {   -logic_name     => 'run_polyphen',
            -module         => 'Bio::EnsEMBL::Variation::Pipeline::ProteinFunction::RunPolyPhen',
            -parameters     => {
                pph_dir => $self->o('pph_dir'),
                @common_params,
            },
            -input_ids      => [],
            -hive_capacity  => 1000,
            -rc_id          => 2,
            -flow_into      => {
                2   => [ 'run_weka' ],
            },
        },
        
        {   -logic_name     => 'run_weka',
            -module         => 'Bio::EnsEMBL::Variation::Pipeline::ProteinFunction::RunWeka',
            -parameters     => { 
                pph_dir => $self->o('pph_dir'),
                @common_params,
            },
            -input_ids      => [],
            -hive_capacity  => 100,
            -rc_id          => 0,
            -flow_into      => {},
        },
        
        {   -logic_name     => 'run_sift',
            -module         => 'Bio::EnsEMBL::Variation::Pipeline::ProteinFunction::RunSift',
            -parameters     => {
                sift_dir    => $self->o('sift_dir'),
                ncbi_dir    => $self->o('ncbi_dir'),
                blastdb     => $self->o('blastdb'),
                @common_params,
            },
            -input_ids      => [],
            -hive_capacity  => 500,
            -rc_id          => 0,
            -flow_into      => {},
        },

        {   -logic_name     => 'rebuild_polyphen_indexes',
            -module         => 'Bio::EnsEMBL::Variation::Pipeline::RebuildIndexes',
            -parameters     => {
                tables  => [ qw(polyphen_prediction polyphen_supplementary_data) ],
                @common_params,
            },
           -input_ids      => [],
            -hive_capacity  => 1,
            -rc_id          => 3,
            -wait_for       => [ 'run_polyphen', 'run_weka' ],
            -flow_into      => {},
        },

        {   -logic_name     => 'rebuild_sift_indexes',
            -module         => 'Bio::EnsEMBL::Variation::Pipeline::RebuildIndexes',
            -parameters     => {
                tables  => [ qw(polyphen_prediction polyphen_supplementary_data) ],
                @common_params,
            },
            -input_ids      => [],
            -hive_capacity  => 1,
            -rc_id          => 3,
            -wait_for       => [ 'run_sift' ],
            -flow_into      => {},
        },

    ];
}

1;

