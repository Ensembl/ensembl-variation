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

=cut
package Bio::EnsEMBL::Variation::Pipeline::ProteinFunction::ProteinFunction_conf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');

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

        debug_mode              => 0,

        # the location of your ensembl checkout, the hive looks here for SQL files etc.

        ensembl_cvs_root_dir    => $ENV{'HOME'} . '/src',
        hive_root_dir           => $ENV{'HOME'} . '/src/ensembl-hive',
        
        pipeline_name           => 'protein_function',
        pipeline_dir            => '/hps/nobackup/production/ensembl/'.$ENV{USER}.'/'.$self->o('pipeline_name'),
        
        species_dir             => $self->o('pipeline_dir').'/'.$self->o('species'),
        
        # directory used for the hive's own output files

        output_dir              => $self->o('species_dir').'/hive_output',

        # this registry file should contain connection details for the core, variation
        # and compara databases (if you are using compara alignments). If you are
        # doing an UPDATE run for either sift or polyphen, the variation database
        # should have existing predictions in the protein_function_predictions table

        ensembl_registry        => $self->o('species_dir').'/ensembl.registry',

        # peptide sequences for all unique translations for this species will be dumped to this file

        fasta_file              => $self->o('species_dir').'/'.$self->o('species').'_translations.fa',
        
        # set this flag to include LRG translations in the analysis

        include_lrg             => 0,

        # include RefSeq transcripts, and edit with accompanying BAM?
        include_refseq          => 0,
        bam                     => '/nfs/production/panda/ensembl/variation/data/dump_vep/interim_GRCh38.p10_knownrefseq_alignments_2017-01-13.bam',

        # GRCh37 bam
        # bam                     => '/nfs/production/panda/ensembl/variation/data/dump_vep/interim_GRCh37.p13_knownrefseq_alignments_2017-01-13.bam',        
        
        # connection details for the hive's own database

        pipeline_db => {
            -host   => 'mysql-ens-var-prod-1.ebi.ac.uk',
            -port   => 4449,
            -user   => 'ensadmin',
            -pass   => $self->o('password'),            
            -dbname => $ENV{USER}.'_'.$self->o('pipeline_name').'_'. $self->o('species') .'_hive',
            -driver => 'mysql',
        },
        
        # configuration for the various resource options used in the pipeline
        
        default_lsf_options => '-qproduction-rh7 -R"select[mem>2000] rusage[mem=2000]" -M2000',
        medmem_lsf_options  => '-qproduction-rh7 -R"select[mem>8000] rusage[mem=8000]" -M8000',
        urgent_lsf_options  => '-qproduction-rh7 -R"select[mem>2000] rusage[mem=2000]" -M2000',
        highmem_lsf_options => '-qproduction-rh7 -R"select[mem>16000] rusage[mem=16000]" -M16000', # this is Sanger LSF speak for "give me 15GB of memory"
        long_lsf_options    => '-qproduction-rh7 -R"select[mem>2000] rusage[mem=2000]" -M2000',

        # Polyphen specific parameters

        # location of the software

        pph_dir                 => '/nfs/production/panda/ensembl/variation/software/polyphen-2.2.2',

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

        sift_dir                => '/nfs/panda/ensemblgenomes/external/sift',

        sift_working            => $self->o('species_dir').'/sift_working',
        
        # the location of blastpgp etc.

        ncbi_dir                => '/nfs/panda/ensemblgenomes/external/ncbi-blast-2+/bin',
        
        # the protein database used to build alignments if you're not using compara

        blastdb                 => '/nfs/production/panda/ensembl/variation/data/sift5.2.2/uniref90/uniref90.fasta',

        # the following parameters mean the same as for polyphen

        sift_run_type           => UPDATE,

        sift_use_compara        => 0,

        sift_max_workers        => 500,
    };
}

#sub pipeline_create_commands {
#    my ($self) = @_;
#    return [
#        'mysql '.$self->dbconn_2_mysql('pipeline_db', 0).q{-e 'DROP DATABASE IF EXISTS }.$self->o('pipeline_db', '-dbname').q{'},
#        @{$self->SUPER::pipeline_create_commands}, 
#        'mysql '.$self->dbconn_2_mysql('pipeline_db', 1).q{-e 'INSERT INTO meta (meta_key, meta_value) VALUES ("hive_output_dir", "}.$self->o('output_dir').q{")'},
#    ];
#}

sub resource_classes {
    my ($self) = @_;
    return {
          'default' => { 'LSF' => $self->o('default_lsf_options') },
          'medmem'  => { 'LSF' => $self->o('medmem_lsf_options')  },
          'urgent'  => { 'LSF' => $self->o('urgent_lsf_options')  },
          'highmem' => { 'LSF' => $self->o('highmem_lsf_options') },
          'long'    => { 'LSF' => $self->o('long_lsf_options')    },
    };
}

sub pipeline_analyses {
    my ($self) = @_;
    
    my @common_params = (
        fasta_file          => $self->o('fasta_file'),
        ensembl_registry    => $self->o('ensembl_registry'),
        species             => $self->o('species'),
        debug_mode          => $self->o('debug_mode'),
    );

    return [
        {   -logic_name => 'init_jobs',
            -module     => 'Bio::EnsEMBL::Variation::Pipeline::ProteinFunction::InitJobs',
            -parameters => {
                sift_run_type   => $self->o('sift_run_type'),
                pph_run_type    => $self->o('pph_run_type'),
                include_lrg     => $self->o('include_lrg'),
                polyphen_dir    => $self->o('pph_dir'),
                sift_dir        => $self->o('sift_dir'),                
                blastdb         => $self->o('blastdb'),
                include_refseq  => $self->o('include_refseq'),
                bam             => $self->o('bam'),
                species_dir     => $self->o('species_dir'),
                @common_params,
            },
            -input_ids  => [{}],
            -rc_name    => 'highmem',
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
            -rc_name        => 'highmem',
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
            -rc_name        => 'default',
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
            -failed_job_tolerance => 10,
            -max_retry_count => 5,
            -input_ids      => [],
            -hive_capacity  => $self->o('sift_max_workers'),
            -rc_name        => 'medmem',
            -flow_into      => {
              -1 => ['run_sift_highmem'],
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

    ];
}

1;

