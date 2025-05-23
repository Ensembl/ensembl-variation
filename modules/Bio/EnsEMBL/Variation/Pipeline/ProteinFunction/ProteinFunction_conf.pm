=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2025] EMBL-European Bioinformatics Institute

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

=cut
package Bio::EnsEMBL::Variation::Pipeline::ProteinFunction::ProteinFunction_conf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');

use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
use Bio::EnsEMBL::Variation::Pipeline::ProteinFunction::Constants qw(FULL UPDATE NONE);

sub default_options {
    my ($self) = @_;

    return {

        # NB: You can find some documentation on this pipeline on confluence here:
        #
        # http://www.ebi.ac.uk/seqdb/confluence/display/EV/Protein+function+pipeline

        # Pipeline wide settings

        # If the debug_mode flag is set to 1 then we will only run the pipeline on a single gene
        # (currently set to BRCA1 in InitJobs.pm), this is useful when testing new installations
        # of sift and polyphen, alterations to the pipeline etc. When set to 0 the full pipeline
        # will be run
        
        hive_force_init => 1,
        hive_use_param_stack => 0,
        hive_use_triggers => 0,
        hive_auto_rebalance_semaphores => 0, 
        hive_no_init => 0,
        hive_default_max_retry_count => 0,
        hive_debug_init => 1,
        debug_mode => 0,

        # the location of your ensembl checkout, the hive looks here for SQL files etc.

        ensembl_cvs_root_dir    => $ENV{'HOME'} . '/src',
        hive_root_dir           => $ENV{'HOME'} . '/src/ensembl-hive',
        
        pipeline_name           => 'protein_function',
        pipeline_dir            => '/hps/nobackup/flicek/ensembl/variation/'.$ENV{USER}.'/'.$self->o('pipeline_name'),
        species_dir             => $self->o('pipeline_dir').'/'.$self->o('species'),
        # directory used for the hive's own output files

        output_dir              => $self->o('species_dir').'/hive_output',

        # this registry file should contain connection details for the core, variation
        # and compara databases (if you are using compara alignments). If you are
        # doing an UPDATE run for either sift or polyphen, the variation database
        # should have existing predictions in the protein_function_predictions table

        ensembl_registry        => $self->o('species_dir').'/ensembl.registry',

        # a file containing history of datachecks ran potentially used to determine
        # if a datacheck can be skipped 

        history_file            => '/nfs/production/flicek/ensembl/production/datachecks/history/vertebrates.json',

        # output dir where datacheck result will be stored

        dc_outdir               => $self->o('pipeline_dir')."/".$self->o('pipeline_name')."_dc_output",

        # if set, fails the datacheck pipeline job if the datacheck fails
        # can be overwritten when running the pipeline

        failures_fatal          => 1,

        # if set, runs the datachecks analysis jobs
        # can be overwritten when running the pipeline

        run_dc                  => 1,

        # the uri of the database server which stores the database of previous release
        # supported format is mysql://[a_user]@[some_host]:[port_number]/[old_dbname|old_release_number]

        old_server_uri          => undef,

        # peptide sequences for all unique translations for this species will be dumped to this file

        fasta_file              => $self->o('species_dir').'/'.$self->o('species').'_translations.fa',
        # set this flag to include LRG translations in the analysis

        include_lrg             => 0,

        # include RefSeq transcripts, and edit with accompanying BAM?
        include_refseq          => 0,
        bam                     => '/nfs/production/flicek/ensembl/variation/data/dump_vep/GCF_000001405.39_GRCh38.p13_knownrefseq_alns.bam',

        # GRCh37 bam
        # bam                     => '/nfs/production/flicek/ensembl/variation/data/dump_vep/interim_GRCh37.p13_knownrefseq_alignments_2017-01-13.bam',
        
        # connection details for the hive's own database
        hive_db_host    => 'mysql-ens-var-prod-2.ebi.ac.uk',
        hive_db_port    => 4521,
        hive_db_user    => 'ensadmin',
        hive_db_name    => $ENV{'USER'} . '_ehive_' . $self->o('pipeline_name') . '_' . $self->o('ensembl_release') . '_' . $self->o('assembly') . '_' . $self->o('species'),

        pipeline_db => {
            -host   => $self->o('hive_db_host'),
            -port   => $self->o('hive_db_port'),
            -user   => $self->o('hive_db_user'),
            -pass   => $self->o('password'),            
            -dbname => $self->o('hive_db_name'),
            -driver => 'mysql',
        },
        
        # configuration for the various resource options used in the pipeline
        
        default_lsf_options => '-qproduction -R"select[mem>2000] rusage[mem=2000]" -M2000',
        medmem_lsf_options  => '-qproduction -R"select[mem>8000] rusage[mem=8000]" -M8000',
        highmem_lsf_options => '-qproduction -R"select[mem>24000] rusage[mem=24000]" -M24000',

        default_slurm_options      => '--partition=production --time=1:30:00 --mem=2G',
        default_long_slurm_options => '--partition=production --time=6:00:00 --mem=2G',
        medmem_slurm_options       => '--partition=production --time=6:00:00 --mem=8G',
        medmem_long_slurm_options  => '--partition=production --time=24:00:00 --mem=8G',
        highmem_slurm_options      => '--partition=production --time=2:30:00 --mem=24G',
        highmem_med_slurm_options  => '--partition=production --time=24:00:00 --mem=24G',
        highmem_long_slurm_options => '--partition=production --time=120:00:00 --mem=24G',

        # Polyphen specific parameters

        # location of the software

        pph_dir                 => '/hps/software/users/ensembl/variation/polyphen-2.2.3',
        pph_data                => $self->o('pph_dir').'/data',

        # where we will keep polyphen's working files etc. as the pipeline runs

        pph_working             => $self->o('species_dir').'/polyphen_working',

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

        pph_run_type            => NONE,

        # set this flag to use compara protein families as the alignments rather than
        # polyphen's own alignment pipeline

        pph_use_compara         => 0,

        # the maximum number of workers to run in parallel for polyphen and weka. Weka 
        # runs much faster then polyphen so you don't need as many workers.

        pph_max_workers         => 500,

        weka_max_workers        => 20,

        # Sift specific parameters
    
        # location of the software

        sift_dir                => '/hps/software/users/ensembl/variation/sift6.2.1',

        sift_working            => $self->o('species_dir').'/sift_working',
        
        # the location of blastpgp etc.

        ncbi_dir                => '/hps/software/users/ensembl/ensw/C8-MAR21-sandybridge/linuxbrew/Cellar/blast/2.2.30/bin',
        
        # the protein database used to build alignments if you're not using compara

        variation_data          => '/nfs/production/flicek/ensembl/variation/data',

        blastdb                 => $self->o('variation_data').'/uniref90/uniref90.fasta',

        # the following parameters mean the same as for polyphen

        sift_run_type           => NONE,

        sift_use_compara        => 0,

        sift_max_workers        => 500,

        dbnsfp_run_type         => NONE,
        dbnsfp_max_workers      => 50,
        dbnsfp_working          => $self->o('species_dir').'/dbnsfp_working',
        dbnsfp_annotation       => { GRCh37 =>
                                      { file => $self->o('variation_data') . '/dbNSFP/4.9c/dbNSFP4.9c_grch37.gz',
                                        version => '4.9c',
                                      },
                                     GRCh38 =>
                                      { file => $self->o('variation_data') . '/dbNSFP/4.9c/dbNSFP4.9c_grch38.gz',
                                        version => '4.9c',
                                      } 
                                    },
        cadd_run_type         => NONE,
        cadd_max_workers      => 50,
        cadd_working          => $self->o('species_dir').'/cadd_working',
        cadd_annotation       => { GRCh37 =>
                                      { file => $self->o('variation_data') . '/CADD/v1.7/grch37/CADD_GRCh37_1.7_whole_genome_SNVs.tsv.gz',
                                        version => 'v1.7',
                                      },
                                   GRCh38 =>
                                      { file => $self->o('variation_data') . '/CADD/v1.7/grch38/CADD_GRCh38_1.7_whole_genome_SNVs.tsv.gz',
                                        version => 'v1.7',
                                      } 
                                  },
    };
}

sub pipeline_create_commands {
  my ($self) = @_;
  return [
    @{$self->SUPER::pipeline_create_commands},  # inheriting database and hive tables' creation
    $self->db_cmd('CREATE TABLE IF NOT EXISTS failure_reason (
        translation_md5 char(32) NOT NULL,
        error_msg varchar(255) NOT NULL,
        analysis char(32) NOT NULL,
        PRIMARY KEY (translation_md5),
        UNIQUE KEY md5_error_analysis  (translation_md5, error_msg, analysis)
        ) ENGINE=InnoDB DEFAULT CHARSET=latin1;'),
  ];
}

sub resource_classes {
    my ($self) = @_;
    return {
          'default'      => { 'LSF' => $self->o('default_lsf_options'),
                              'SLURM' => $self->o('default_slurm_options') },
          'default_long' => { 'LSF' => $self->o('default_lsf_options'),
                              'SLURM' => $self->o('default_long_slurm_options') },
          'medmem'       => { 'LSF' => $self->o('medmem_lsf_options'),
                              'SLURM' => $self->o('medmem_slurm_options')  },
          'medmem_long'  => { 'LSF' => $self->o('medmem_lsf_options'),
                              'SLURM' => $self->o('medmem_long_slurm_options')  },
          'highmem'      => { 'LSF' => $self->o('highmem_lsf_options'),
                              'SLURM' => $self->o('highmem_slurm_options') },
          'highmem_med'  => { 'LSF' => $self->o('highmem_lsf_options'),
                              'SLURM' => $self->o('highmem_med_slurm_options') },
          'highmem_long' => { 'LSF' => $self->o('highmem_lsf_options'),
                              'SLURM' => $self->o('highmem_long_slurm_options') },
    };
}

sub pipeline_analyses {
    my ($self) = @_;
    
    my @common_params = (
        fasta_file          => $self->o('fasta_file'),
        ensembl_registry    => $self->o('ensembl_registry'),
        species             => $self->o('species'),
        ensembl_release     => $self->o('ensembl_release'),
        assembly            => $self->o('assembly'),
        debug_mode          => $self->o('debug_mode'),
    );

    return [
        {   -logic_name => 'init_jobs',
            -module     => 'Bio::EnsEMBL::Variation::Pipeline::ProteinFunction::InitJobs',
            -parameters => {
                sift_run_type   => $self->o('sift_run_type'),
                pph_run_type    => $self->o('pph_run_type'),
                dbnsfp_run_type => $self->o('dbnsfp_run_type'),
                dbnsfp_working  => $self->o('dbnsfp_working'),
                dbnsfp_annotation => $self->o('dbnsfp_annotation'),
                cadd_run_type   => $self->o('cadd_run_type'),
                cadd_working    => $self->o('cadd_working'),
                cadd_annotation => $self->o('cadd_annotation'),
                include_lrg     => $self->o('include_lrg'),
                polyphen_dir    => $self->o('pph_dir'),
                polyphen_data   => $self->o('pph_data'),
                sift_dir        => $self->o('sift_dir'),                
                blastdb         => $self->o('blastdb'),
                include_refseq  => $self->o('include_refseq'),
                bam             => $self->o('bam'),
                species_dir     => $self->o('species_dir'),
                use_compara     => $self->o('sift_use_compara'),
                run_dc          => $self->o('run_dc'),
                old_server_uri  => $self->o('old_server_uri'),
                @common_params,
            },
            -input_ids  => [{}],
            -rc_name    => 'highmem',
            -max_retry_count => 0,
            -flow_into  => {
                '2->A' => [ 'run_polyphen' ],
                '3->A' => [ 'run_sift' ],
                '4->A' => [ 'run_dbnsfp' ],
                '5->A' => [ 'run_cadd' ],
                'A->1' => WHEN(
                    '#run_dc#' => [ 'datacheck' ]
                )
            },
        },

        {   -logic_name     => 'run_polyphen',
            -module         => 'Bio::EnsEMBL::Variation::Pipeline::ProteinFunction::RunPolyPhen',
            -parameters     => {
                pph_dir     => $self->o('pph_dir'),
                pph_data    => $self->o('pph_data'),
                pph_working => $self->o('pph_working'),
                use_compara => $self->o('pph_use_compara'),
                @common_params,
            },
            -max_retry_count => 0,
            -input_ids      => [],
            -hive_capacity  => $self->o('pph_max_workers'),
            -rc_name        => 'highmem_med',
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
                pph_data        => $self->o('pph_data'),
                @common_params,
            },
            -max_retry_count => 0,
            -input_ids      => [],
            -hive_capacity  => $self->o('weka_max_workers'),
            -rc_name        => 'default',
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
            -failed_job_tolerance => 10,
            -max_retry_count => 0,
            -input_ids      => [],
            -hive_capacity  => $self->o('sift_max_workers'),
            -rc_name        => 'medmem_long',
            -flow_into      => {
                -1 => ['run_sift_highmem']
            }
        },

        {   -logic_name     => 'run_sift_highmem',
            -module         => 'Bio::EnsEMBL::Variation::Pipeline::ProteinFunction::RunSift',
            -parameters     => {
                sift_dir        => $self->o('sift_dir'),
                sift_working    => $self->o('sift_working'),
                ncbi_dir        => $self->o('ncbi_dir'),
                blastdb         => $self->o('blastdb'),
                use_compara     => $self->o('sift_use_compara'),
                @common_params,
            },
            -input_ids      => [],
            -rc_name        => 'highmem',
        },

        {   -logic_name     => 'run_dbnsfp',
            -module         => 'Bio::EnsEMBL::Variation::Pipeline::ProteinFunction::RunDbNSFP',
            -parameters     => {
                dbnsfp_working    => $self->o('dbnsfp_working'),
                dbnsfp_annotation => $self->o('dbnsfp_annotation'),
                @common_params,
            },
            -failed_job_tolerance => 0,
            -max_retry_count => 0,
            -input_ids      => [],
            -hive_capacity  => $self->o('dbnsfp_max_workers'),
            -rc_name        => 'medmem',
        },

        {   -logic_name     => 'run_cadd',
            -module         => 'Bio::EnsEMBL::Variation::Pipeline::ProteinFunction::RunCADD',
            -parameters     => {
                cadd_working    => $self->o('cadd_working'),
                cadd_annotation => $self->o('cadd_annotation'),
                @common_params,
            },
            -failed_job_tolerance => 0,
            -max_retry_count => 0,
            -input_ids      => [],
            -hive_capacity  => $self->o('cadd_max_workers'),
            -rc_name        => 'medmem',
        },

        {   -logic_name      => 'datacheck',
            -module          => 'Bio::EnsEMBL::DataCheck::Pipeline::RunDataChecks',
            -parameters      => {
                datacheck_names => [
                    'CompareProteinFunctionPredictions',
                    'ProteinFunctionPredictions'
                ],
                registry_file  => $self->o('ensembl_registry'),
                history_file   => $self->o('history_file'),
                output_dir     => $self->o('dc_outdir'),
                failures_fatal => $self->o('failures_fatal'),
                @common_params
            },            
            -input_ids            => [], #default
            -hive_capacity        => 1,
            -analysis_capacity    => 1,
            -rc_name              => 'default_long',
            -failed_job_tolerance => 0,
            -max_retry_count      => 0,
        },

    ];
}

1;
