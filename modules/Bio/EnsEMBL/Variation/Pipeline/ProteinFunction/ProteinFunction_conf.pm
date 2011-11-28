package Bio::EnsEMBL::Variation::Pipeline::ProteinFunction::ProteinFunction_conf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');

use Bio::EnsEMBL::Variation::Pipeline::ProteinFunction::Constants qw(FULL UPDATE NONE);

sub default_options {
    my ($self) = @_;

    return {
    
        ensembl_cvs_root_dir    => $ENV{'HOME'}.'/ensembl-branches/HEAD/',

        pipeline_name           => 'protein_function',

        pipeline_dir            => '/lustre/scratch101/ensembl/gr5/protein_function',
        
        output_dir              => $self->o('pipeline_dir').'/hive_output',

        ensembl_registry        => $self->o('pipeline_dir').'/ensembl.registry',

        proteins_fasta          => $self->o('pipeline_dir').'/ensembl_65_proteins.fa',

        species                 => 'human',

        pph_dir                 => $self->o('pipeline_dir').'/polyphen-2.1.0',
        
        sift_dir                => $self->o('pipeline_dir').'/sift-4.0.3',
        
        ncbi_dir                => '/software/ncbiblast/bin',
        
        blastdb                 => $self->o('pipeline_dir').'/blastdb/swissprot_trembl.uni',

        sift_run_type           => FULL,

        polyphen_run_type       => FULL,

        sift_use_compara        => 1,
        
        pph_use_compara         => 1,

        include_lrg             => 0,
        
        use_existing_table      => 0,

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
        2 => { -desc => 'highmem',  'LSF' => '-R"select[mem>8000] rusage[mem=8000]" -M8000000 -q long'},
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
                sift_run_type       => $self->o('sift_run_type'),
                polyphen_run_type   => $self->o('polyphen_run_type'),
                include_lrg         => $self->o('include_lrg'),
                use_existing_table  => $self->o('use_existing_table'),
                @common_params,
            },
            -input_ids  => [{}],
            -rc_id      => 2,
            -flow_into  => {
                2 => [ 'run_polyphen' ],
                3 => [ 'run_sift' ],
            },
        },

        {   -logic_name     => 'run_polyphen',
            -module         => 'Bio::EnsEMBL::Variation::Pipeline::ProteinFunction::RunPolyPhen',
            -parameters     => {
                pph_dir     => $self->o('pph_dir'),
                use_compara => $self->o('pph_use_compara'),
                @common_params,
            },
            -max_retry_count => 0,
            -input_ids      => [],
            -hive_capacity  => 500,
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
            -max_retry_count => 0,
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
                use_compara => $self->o('sift_use_compara'),
                @common_params,
            },
            -max_retry_count => 0,
            -input_ids      => [],
            -hive_capacity  => 500,
            -rc_id          => 0,
            -flow_into      => {},
        },
    ];
}

1;

