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
    
        # the location of your ensembl checkout, the hive looks here for SQL files etc.

        ensembl_cvs_root_dir    => $ENV{'HOME'}.'/ensembl-branches/HEAD/',

        pipeline_name           => 'protein_function',

        pipeline_dir            => '/lustre/scratch101/ensembl/'.$ENV{USER}.'/'.$self->o('pipeline_name'),
        
        # directory used for the hive's own output files

        output_dir              => $self->o('pipeline_dir').'/hive_output',

        # this registry file should contain conection details for the core, variation
        # and compara databases (if you are using compara alignments). If you are
        # doing an UPDATE run for either sift or polyphen, the variation database
        # should have existing predictions in the protein_function_predictions table

        ensembl_registry        => $self->o('pipeline_dir').'/ensembl.registry',

        # peptide sequences for all unique translations will be dumped to this file

        fasta_file              => $self->o('pipeline_dir').'/required_translations.fa',
        
        # set this flag to include LRG translations in the analysis

        include_lrg             => 1,
        
        # connection details for the hive's own database

        pipeline_db => {
            -host   => 'ens-variation',
            -port   => 3306,
            -user   => 'ensadmin',
            -pass   => $self->o('password'),            
            -dbname => $ENV{USER}.'_'.$self->o('pipeline_name').'_hive',
        },
        
        # configuration for the various resource options used in the pipeline
        
        default_lsf_options => '-R"select[mem>2000] rusage[mem=2000]" -M2000000',
        urgent_lsf_options  => '-q yesterday -R"select[mem>2000] rusage[mem=2000]" -M2000000',
        highmem_lsf_options => '-q long -R"select[mem>8000] rusage[mem=8000]" -M8000000',
        medmem_lsf_options  => '-R"select[mem>4000] rusage[mem=4000]" -M4000000',
        long_lsf_options    => '-q long -R"select[mem>2000] rusage[mem=2000]" -M2000000',

        # Polyphen specific parameters

        # location of the software

        pph_dir                 => '/software/ensembl/variation/polyphen-2.2.2',

        # where we will keep polyphen's working files etc. as the pipeline runs

        pph_working             => $self->o('pipeline_dir').'/polyphen_working',
        
        # specify the Weka classifier models here, if you don't want predictions from 
        # one of the classifier models set the value to the empty string

        humdiv_model            => $self->o('pph_dir').'/models/HumDiv.UniRef100.NBd.f11.model',
        
        humvar_model            => $self->o('pph_dir').'/models/HumVar.UniRef100.NBd.f11.model',

        # the run type for polyphen (& sift) can be one of FULL to run predictions for
        # all translations regardless of whether we already have predictions in the
        # database, NONE to exclude this analysis, or UPDATE to run predictions for any
        # new or changed translations in the database. The variation database specified 
        # in the registry above is used to identify translations we already have 
        # predictions for.

        pph_run_type            => FULL,

        # set this flag to use compara protein families as the alignments rather than
        # polyphen's own alignment pipeline

        pph_use_compara         => 0,

        # the maximum number of workers to run in parallel for polyphen and weka. Weka 
        # runs much faster then polyphen so you don't need as many workers.

        pph_max_workers         => 500,

        weka_max_workers        => 20,

        # Sift specific parameters
    
        # location of the software

        sift_dir                => '/software/ensembl/variation/sift4.0.5',

        sift_working            => $self->o('pipeline_dir').'/sift_working',
        
        # the location of blastpgp etc.

        ncbi_dir                => '/software/ncbiblast/bin',
        
        # the protein database used to build alignments if you're not using compara

        blastdb                 => '/data/blastdb/Ensembl/variation/sift4.0.5/uniprot/swiss_trembl.uni',

        # the following parameters mean the same as for polyphen

        sift_run_type           => FULL,

        sift_use_compara        => 0,

        sift_max_workers        => 500,
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
        4 => { -desc => 'medmem',   'LSF' => $self->o('medmem_lsf_options')    },
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
            -hive_capacity  => $self->o('pph_max_workers'),
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
            -hive_capacity  => $self->o('weka_max_workers'),
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
            -hive_capacity  => $self->o('sift_max_workers'),
            -rc_id          => 4,
            -flow_into      => {},
        },
    ];
}

1;

