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

package Bio::EnsEMBL::Variation::Pipeline::VariationConsequence_full_conf;

use strict;
use warnings;
use File::Spec::Functions qw(catfile catdir);

use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf qw(WHEN ELSE);
use Bio::EnsEMBL::Variation::Pipeline::ProteinFunction::Constants qw(FULL UPDATE NONE);

sub default_options {
    my ($self) = @_;

    # The hash returned from this function is used to configure the
    # pipeline, you can supply any of these options on the command
    # line to override these default values.
    
    # You shouldn't need to edit anything in this file other than
    # these values, if you find you do need to then we should probably
    # make it an option here, contact the variation team to discuss
    # this - patches are welcome!

    return {
      
        # List of species to use
        species => [],

        # general pipeline options that you should change to suit your environment
       
        hive_force_init => 1,
        hive_use_param_stack => 0,
        hive_use_triggers => 0,
        hive_auto_rebalance_semaphores => 0, 
        hive_no_init => 0,
        # the location of your checkout of the ensembl API (the hive looks for SQL files here)
        
        debug_mode              => 0,
        
        ensembl_cvs_root_dir    => $ENV{'HOME'} . '/src',
        hive_root_dir           => $self->o('ensembl_cvs_root_dir') . '/ensembl-hive', 
        # a name for your pipeline (will also be used in the name of the hive database)
        
        pipeline_name           => 'variation_consequence_all',

        # a directory to keep hive output files and your registry file, you should
        # create this if it doesn't exist

        pipeline_dir            => '/hps/nobackup/production/ensembl/' . $ENV{USER} . '/' . $self->o('pipeline_name'),
        
        species_dir             => $self->o('pipeline_dir').'/#species#',
        output_dir              => $self->o('species_dir').'/hive_output',

        # a standard ensembl registry file containing connection parameters
        # for your target database(s) (and also possibly aliases for your species
        # of interest that you can then supply to init_pipeline.pl with the -species
        # option)
        
        reg_file                => $self->o('pipeline_dir').'#species#/ensembl.registry',

        # points to a FASTA file, much faster than using DB for sequence lookup if available
        fasta_dir => undef,
        fasta_file => $self->o('fasta_dir') ? catfile($self->o('fasta_dir'), '#species#', '#species#.fa') : undef,

        # if set to 1 this option tells the transcript_effect analysis to disambiguate
        # ambiguity codes in single nucleotide alleles, so e.g. an allele string like
        # 'T/M' will be treated as if it were 'T/A/C' (this was a request from ensembl
        # genomes and we don't use it by default in the ensembl variation pipeline)
        
        disambiguate_single_nucleotide_alleles => 0,

        # configuration for the various resource options used in the pipeline
        # Users of other farms should change these here, or override them on
        # the command line to suit your farm. The names of each option hopefully
        # reflect their usage, but you may want to change the details (memory
        # requirements, queue parameters etc.) to suit your own data
        
        default_lsf_options => '-qproduction-rh7 -R"select[mem>2000] rusage[mem=2000]" -M2000',
        medmem_lsf_options  => '-qproduction-rh7 -R"select[mem>4000] rusage[mem=4000]" -M4000',
        urgent_lsf_options  => '-qproduction-rh7 -R"select[mem>2000] rusage[mem=2000]" -M2000',
        highmem_lsf_options => '-qproduction-rh7 -R"select[mem>15000] rusage[mem=15000] span[hosts=1]" -M15000 -n4', # this is LSF speak for "give me 15GB of memory"
        long_lsf_options    => '-qproduction-rh7 -R"select[mem>2000] rusage[mem=2000]" -M2000',

        # options controlling the number of workers used for the parallelisable analyses
        # these default values seem to work for most species

        transcript_effect_capacity      => 50,
        set_variation_class_capacity    => 10,
        
        # set this flag to 1 to include LRG transcripts in the transcript effect analysis

        include_lrg => 1, 

        # set this flag to 1 to try and identify genetic markers in SetVariationClass module
        # This is very specific to data imported from dbSNP by ensembl.
        # ensembl genomes might need different methods for idenfifying markers:
        # for the future add identify_marker_eg flag and add code to SetVariationClass module 
        identify_marker_e => 1, 

        # Limit analysis to specific gene biotypes
        limit_biotypes => [],

        # create MTMP_transcript_variation
        mtmp_table => 1,

        # sort variation_feature before we start?
        # disable this if you are sure the table is already sorted
        # or if the table is sufficiently small that it won't make much difference
        sort_variation_feature => 0,

        # sets the maximum distance to a transcript for which up/downstream consequences are assessed
        max_distance => undef,

        # these flags control which parts of the pipeline are run

        run_transcript_effect   => 1,
        run_variation_class     => 1,

        # connection parameters for the hive database, you should supply the hive_db_password
        # option on the command line to init_pipeline.pl (parameters for the target database
        # should be set in the registry file defined above)

        # Should hive use triggeres?
        hive_use_triggers       => 0,

        # init_pipeline.pl will create the hive database on this machine, naming it
        # <username>_<pipeline_name>, and will drop any existing database with this
        # name

        hive_db_host    => 'mysql-ens-var-prod-1',
        hive_db_port    => 4449,
        hive_db_user    => 'ensadmin',

        pipeline_db => {
            -host   => $self->o('hive_db_host'),
            -port   => $self->o('hive_db_port'),
            -user   => $self->o('hive_db_user'),
            -pass   => $self->o('hive_db_password'),            
            -dbname => $ENV{'USER'}.'_'.$self->o('pipeline_name'),
            -driver => 'mysql',
        },
        
        # PROTEIN FEATURES PART
        # set this flag to include LRG translations in the analysis
        include_lrg             => 0,

        # include RefSeq transcripts, and edit with accompanying BAM?
        include_refseq          => 0,
        bam                     => '/nfs/production/panda/ensembl/variation/data/dump_vep/interim_GRCh38.p10_knownrefseq_alignments_2017-01-13.bam',
        
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

sub resource_classes {
    my ($self) = @_;
    return {
          'default' => { 'LSF' => $self->o('default_lsf_options') },
          'urgent'  => { 'LSF' => $self->o('urgent_lsf_options')  },
          'highmem' => { 'LSF' => $self->o('highmem_lsf_options') },
          'long'    => { 'LSF' => $self->o('long_lsf_options')    },
          'medmem'  => { 'LSF' => $self->o('medmem_lsf_options') },
    };
}

sub beekeeper_extra_cmdline_options {
    my $self = shift;
    return "-reg_conf " . $self->o("reg_file");
}

sub pipeline_wide_parameters {
  my ($self) = @_;

  return {
    %{$self->SUPER::pipeline_wide_parameters},
        run_transcript_effect   => $self->o('run_transcript_effect'),
        run_variation_class     => $self->o('run_variation_class'),
  };
}
                        

sub pipeline_analyses {
  my ($self) = @_;

  my @common_params = (
    ensembl_registry    => $self->o('reg_file'),
    pipeline_dir => catdir($self->o('pipeline_dir'), '#species#'),
    fasta_file          => $self->o('fasta_file'),
    debug_mode          => $self->o('debug_mode'),
  );

  my @rebuild_tables = qw(transcript_variation variation_hgvs variation_genename);
  push @rebuild_tables, 'MTMP_transcript_variation' if $self->o('mtmp_table');

  my @analyses;
  push @analyses, (
        {   -logic_name => 'species_factory',
            -module     => 'Bio::EnsEMBL::Production::Pipeline::Production::SpeciesFactory',
            -parameters => {
                db_types => [ 'variation' ],
                species  => $self->o('species'),
            },
            -meadow_type       => 'LOCAL',
            -input_ids  => [{}],
            -rc_name    => 'default',
            -flow_into  => {
                2 => [ 'init_jobs' ],
            },
        },
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
            -rc_name    => 'highmem',
            -flow_into  => {
                '2->A' => [ 'run_polyphen' ],
                '3->A' => [ 'run_sift' ],
                'A->1' => [ 'protein_function_cleanup' ],
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
            -max_retry_count => 0,
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

        {   -logic_name => 'protein_function_cleanup',
            -module     => 'Bio::EnsEMBL::Variation::Pipeline::ProteinFunction::Cleanup',
            -parameters => {
                @common_params,
            },
            -meadow_type       => 'LOCAL',
            -max_retry_count => 0,
            -rc_name    => 'default',
            -flow_into      => {
              1 => ['transcript_effect_start'],
            }
        },
      );
      
  push @analyses, (
    {   -logic_name => 'transcript_effect_start',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -meadow_type       => 'LOCAL',
      -max_retry_count => 0,
      -rc_name    => 'default',
      -flow_into  => {
        1 => WHEN('#run_transcript_effect#' => 'init_transcript_effect',
             ELSE 'init_variation_class'),
      },
    },
  );

  push @analyses, (
    {   -logic_name => 'init_transcript_effect',
      -module     => 'Bio::EnsEMBL::Variation::Pipeline::InitTranscriptEffect',
      -parameters => {
        include_lrg => $self->o('include_lrg'),
        limit_biotypes => $self->o('limit_biotypes'),
        mtmp_table => $self->o('mtmp_table'),
        sort_variation_feature => $self->o('sort_variation_feature'),
        @common_params,
      },
      -hive_capacity  => 5,
      -max_retry_count => 0,
      -input_ids  => [],
      -rc_name    => 'long',
      -flow_into  => {
        '2->A' => [ 'transcript_effect' ],
        '3->B' => [ 'transcript_effect_highmem' ],
        'A->1' => [ 'finish_transcript_effect' ],
        'B->1' => [ 'finish_transcript_effect' ],
      },
    },
    {   -logic_name     => 'transcript_effect',
      -module         => 'Bio::EnsEMBL::Variation::Pipeline::TranscriptEffect',
      -parameters     => { 
        disambiguate_single_nucleotide_alleles => $self->o('disambiguate_single_nucleotide_alleles'),
        mtmp_table => $self->o('mtmp_table'),
        max_distance => $self->o('max_distance'),
        @common_params,
      },
      -max_retry_count => 0,
      -input_ids      => [],
      -hive_capacity  => $self->o('transcript_effect_capacity'),
      -rc_name        => 'default',
      -flow_into      => {
        -1 => ['transcript_effect_highmem'],
      }
    },

    {   -logic_name     => 'transcript_effect_highmem',
      -module         => 'Bio::EnsEMBL::Variation::Pipeline::TranscriptEffect',
      -parameters     => {
        disambiguate_single_nucleotide_alleles => $self->o('disambiguate_single_nucleotide_alleles'),
        mtmp_table => $self->o('mtmp_table'),
        max_distance => $self->o('max_distance'),
        @common_params,
      },
      -max_retry_count => 0,
      -input_ids      => [],
      -hive_capacity  => $self->o('transcript_effect_capacity'),
      -rc_name        => 'highmem',
      -can_be_empty   => 1,
    },

    {   -logic_name     => 'finish_transcript_effect',
      -module         => 'Bio::EnsEMBL::Variation::Pipeline::FinishTranscriptEffect',
      -parameters     => {
        @common_params,
      },
      -max_retry_count => 0,
      -input_ids      => [],
      -hive_capacity  => 5,
      -rc_name        => 'highmem',
      -flow_into      => {
        1 => ['rebuild_tv_indexes'],
      },
    },

    {   -logic_name     => 'rebuild_tv_indexes',
      -module         => 'Bio::EnsEMBL::Variation::Pipeline::RebuildIndexes',
      -parameters     => {
        tables => \@rebuild_tables,
        @common_params,
      },
      -max_retry_count => 0,
      -input_ids      => [],
      -hive_capacity  => 5,
      -rc_name        => 'urgent',
      -flow_into      => {
        1 => ['check_transcript_variation', 'update_variation_feature'],
      },
    },

    {   -logic_name     => 'check_transcript_variation',
      -module         => 'Bio::EnsEMBL::Variation::Pipeline::CheckTranscriptVariation',
      -parameters     => {
        @common_params,
      },
      -max_retry_count => 0,
      -input_ids      => [],
      -hive_capacity  => 5,
      -rc_name        => 'default',
      -flow_into      => {},
    },

    {   -logic_name     => 'update_variation_feature',
      -module         => 'Bio::EnsEMBL::Variation::Pipeline::UpdateVariationFeature',
      -parameters     => {
        @common_params,
      },
      -max_retry_count => 0,
      -input_ids      => [],
      -hive_capacity  => 5,
      -rc_name        => 'urgent',
      -flow_into      => {
        1 => WHEN('#run_variation_class#' => 'init_variation_class'),
      },
    }, 
  );

  push @analyses, (
    {   -logic_name     => 'init_variation_class',
      -module         => 'Bio::EnsEMBL::Variation::Pipeline::InitVariationClass',
      -parameters     => {
        num_chunks  => 50,

        @common_params,
      },
      -max_retry_count => 0,
      -input_ids      => [],
      -hive_capacity  => 5,
      -rc_name        => 'default',
      -flow_into      => {
        '1->A' => [ 'set_variation_class' ],
        'A->2' => [ 'finish_variation_class' ],
      },
    },

    {   -logic_name     => 'set_variation_class',
      -module         => 'Bio::EnsEMBL::Variation::Pipeline::SetVariationClass',
      -parameters     => {
        identify_marker_e => $self->o('identify_marker_e'), 
        @common_params,
      },
      -max_retry_count => 0,
      -input_ids      => [],
      -hive_capacity  => $self->o('set_variation_class_capacity'),
      -rc_name        => 'default',
      -flow_into      => {},
    },

    {   -logic_name     => 'finish_variation_class',
      -module         => 'Bio::EnsEMBL::Variation::Pipeline::FinishVariationClass',
      -parameters     => {
        @common_params,
      },
      -max_retry_count => 0,
      -input_ids      => [],
      -hive_capacity  => 5,
      -rc_name        => 'urgent',
      -flow_into      => {},
    },

  );

  return \@analyses;
}

1;

