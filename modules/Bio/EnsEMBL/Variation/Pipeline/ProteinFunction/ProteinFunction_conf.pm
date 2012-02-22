package Bio::EnsEMBL::Variation::Pipeline::ProteinFunction::ProteinFunction_conf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');

use Bio::EnsEMBL::Variation::Pipeline::ProteinFunction::Constants qw(FULL UPDATE NONE);

sub default_options {
    my ($self) = @_;

    return {

        # Pipeline wide commands

        species                 => 'Homo_sapiens',
    
        ensembl_cvs_root_dir    => $ENV{'HOME'}.'/ensembl-branches/HEAD/',

        pipeline_name           => 'protein_function',

        pipeline_dir            => '/lustre/scratch101/ensembl/'.$ENV{USER}.'/protein_function',
        
        output_dir              => $self->o('pipeline_dir').'/hive_output',

        ensembl_registry        => $self->o('pipeline_dir').'/ensembl.registry',

        fasta_file              => $self->o('pipeline_dir').'/required_translations.fa',
        
        include_lrg             => 1,
        
        pipeline_db => {
            -host   => 'ens-variation',
            -port   => 3306,
            -user   => 'ensadmin',
            -pass   => $self->o('password'),            
            -dbname => $ENV{USER}.'_'.$self->o('pipeline_name').'_hive',
        },

        # Polyphen parameters

        pph_dir                 => '/software/ensembl/variation/polyphen-2.2.2',

        pph_working             => $self->o('pipeline_dir').'/polyphen_working',
        
        # if you don't want predictions from one of the classifier models set the value
        # here to the empty string

        humdiv_model            => $self->o('pph_dir').'/models/HumDiv.UniRef100.NBd.f11.model',
        
        humvar_model            => $self->o('pph_dir').'/models/HumVar.UniRef100.NBd.f11.model',

        pph_run_type            => FULL,

        pph_use_compara         => 0,
        
        # Sift parameters

        sift_dir                => '/software/ensembl/variation/sift4.0.5',

        sift_working            => $self->o('pipeline_dir').'/sift_working',
        
        ncbi_dir                => '/software/ncbiblast/bin',
        
        blastdb                 => '/data/blastdb/Ensembl/variation/sift4.0.5/uniprot/swiss_trembl.uni',

        sift_run_type           => FULL,

        sift_use_compara        => 0,
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
        0 => { -desc => 'default',  'LSF' => '-R"select[mem>2000] rusage[mem=2000]" -M2000000' },
        1 => { -desc => 'urgent',   'LSF' => '-q yesterday -R"select[mem>2000] rusage[mem=2000]" -M2000000' },
        2 => { -desc => 'highmem',  'LSF' => '-R"select[mem>8000] rusage[mem=8000]" -M8000000 -q long'},
        3 => { -desc => 'long',     'LSF' => '-q long -R"select[mem>2000] rusage[mem=2000]" -M2000000 -q long' },
        4 => { -desc => 'medmem',   'LSF' => '-R"select[mem>4000] rusage[mem=4000]" -M4000000' },
    };
}

sub pipeline_analyses {
    my ($self) = @_;
    
    my @common_params = (
        fasta_file          => $self->o('fasta_file'),
        ensembl_registry    => $self->o('ensembl_registry'),
        species             => $self->o('species'),
    );

    return [
        {   -logic_name => 'init_jobs',
            -module     => 'Bio::EnsEMBL::Variation::Pipeline::ProteinFunction::InitJobs',
            -parameters => {
                sift_run_type   => $self->o('sift_run_type'),
                pph_run_type    => $self->o('pph_run_type'),
                include_lrg     => $self->o('include_lrg'),
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
                pph_working => $self->o('pph_working'),
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
                pph_dir         => $self->o('pph_dir'),
                humdiv_model    => $self->o('humdiv_model'),
                humvar_model    => $self->o('humvar_model'),
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
                sift_dir        => $self->o('sift_dir'),
                sift_working    => $self->o('sift_working'),
                ncbi_dir        => $self->o('ncbi_dir'),
                blastdb         => $self->o('blastdb'),
                use_compara     => $self->o('sift_use_compara'),
                @common_params,
            },
            -max_retry_count => 0,
            -input_ids      => [],
            -hive_capacity  => 500,
            -rc_id          => 4,
            -flow_into      => {},
        },
    ];
}

1;

