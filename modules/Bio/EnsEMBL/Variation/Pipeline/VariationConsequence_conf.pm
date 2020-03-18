=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2020] EMBL-European Bioinformatics Institute

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

package Bio::EnsEMBL::Variation::Pipeline::VariationConsequence_conf;

use strict;
use warnings;
use Bio::EnsEMBL::Hive::Version 2.5;
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
use base ('Bio::EnsEMBL::Hive::PipeConfig::EnsemblGeneric_conf');

sub default_options {
    my ($self) = @_;

    # The hash returned from this function is used to configure the
    # pipeline, you can supply any of these options on the command
    # line to override these default values.
    
    # You shouldn't need to edit anything in this file other than
    # these values, if you find you do need to then we should probably
    # make it an option here, contact the variation team to discuss
    # this - patches are welcome!
    
    my $login = `whoami`;
    chomp $login;

    return {

        # general pipeline options that you should change to suit your environment
        hive_debug_init => 1, # If set to 1, will show the objects (analyses, data-flow rules, etc) that are parsed from the PipeConfig file.
        hive_default_max_retry_count => 0,
        hive_force_init => 1,
        hive_use_param_stack => 0,
        hive_use_triggers => 0,
        hive_auto_rebalance_semaphores => 0, 
        hive_no_init => 0,
        # the location of your checkout of the ensembl API (the hive looks for SQL files here)
        
        ensembl_cvs_root_dir    => $ENV{'HOME'} . '/src',
        hive_root_dir           => $ENV{'HOME'} . '/src/ensembl-hive', 
        # a name for your pipeline (will also be used in the name of the hive database)
        
        pipeline_name           => 'variation_consequence',

        # a directory to keep hive output files and your registry file, you should
        # create this if it doesn't exist

        pipeline_dir            => '/hps/nobackup/production/ensembl/' . $login . '/' . $self->o('pipeline_name') . '/' . $self->o('species'),

        # a directory where hive workers will dump STDOUT and STDERR for their jobs
        # if you use lots of workers this directory can get quite big, so it's
        # a good idea to keep it on lustre, or some other place where you have a 
        # healthy quota!
        
        output_dir              => $self->o('pipeline_dir').'/hive_output',

        # a standard ensembl registry file containing connection parameters
        # for your target database(s) (and also possibly aliases for your species
        # of interest that you can then supply to init_pipeline.pl with the -species
        # option)
        
        reg_file                => $self->o('pipeline_dir').'/ensembl.registry',

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
        
        default_lsf_options => '-qproduction-rh74 -R"select[mem>2000] rusage[mem=2000]" -M2000',
        medmem_lsf_options  => '-qproduction-rh74 -R"select[mem>5000] rusage[mem=5000]" -M5000',
        highmem_lsf_options => '-qproduction-rh74 -R"select[mem>15000] rusage[mem=15000] span[hosts=1]" -M15000 -n4', # this is LSF speak for "give me 15GB of memory"

        # options controlling the number of workers used for the parallelisable analyses
        # these default values seem to work for most species

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
        sort_variation_feature => 1,

        # points to a FASTA file, much faster than using DB for sequence lookup if available
        fasta => undef,

        # sets the maximum distance to a transcript for which up/downstream consequences are assessed
        max_distance => undef,

        # these flags control which parts of the pipeline are run
        run_transcript_effect   => 1,
        run_variation_class     => 1,

        # Human runs switch off run_var_class and set max_distance to 0 by default. To override
        # this behaviour, set this flag to 1
        human_default_override		=> 0,

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
            -dbname => $ENV{'USER'}.'_'.$self->o('pipeline_name').'_'.$self->o('species'),
            -driver => 'mysql',
            -reconnect_when_lost => 1
        },
    };
}


sub resource_classes {
    my ($self) = @_;
    return {
          'default' => { 'LSF' => $self->o('default_lsf_options') },
          'highmem' => { 'LSF' => $self->o('highmem_lsf_options') },
          'medmem'  => { 'LSF' => $self->o('medmem_lsf_options') },
    };
}

sub pipeline_analyses {
    my ($self) = @_;

    my @common_params = (
        ensembl_registry    => $self->o('reg_file'),
        species             => $self->o('species'),
        pipeline_dir => $self->o('pipeline_dir'),
        max_distance => ($self->o('species') =~ /homo_sapiens|human/ && (! $self->o('human_default_override'))) ? 0 : $self->o('max_distance'),
    );
   
    my @analyses;

    if ($self->o('run_transcript_effect')) {
        push @analyses, (
          { -logic_name => 'init_transcript_effect',
            -module => 'Bio::EnsEMBL::Variation::Pipeline::VariationConsequencePreTasks',
            -parameters => {
              mtmp_table => $self->o('mtmp_table'),
              fasta => $self->o('fasta'),
              sort_variation_feature => $self->o('sort_variation_feature'),
              @common_params,
            },
            -rc_name   => 'default',
            -flow_into => {
              1 => ['gene_factory'],
              2 => ['rebuild_tv_indexes'],
            },
            -input_ids  => [{}],
          },
          { -logic_name => 'gene_factory',
            -module => 'Bio::EnsEMBL::Variation::Pipeline::GeneFactory',
            -hive_capacity  => 50,

            -parameters => {
              include_lrg => $self->o('include_lrg'),
              limit_biotypes => $self->o('limit_biotypes'),
              @common_params,
            },
            -rc_name   => 'default',
            -flow_into => { 
              '2->A' => ['dump_variation_gene_name'], 
              'A->1' => ['web_index_load'], 
            },
          },
          { -logic_name => 'dump_variation_gene_name',
            -module => 'Bio::EnsEMBL::Variation::Pipeline::DumpVariationGeneName',
            -hive_capacity  => 50,
            -parameters => {
              @common_params,
            },
            -rc_name   => 'default',
            -flow_into => { 
              -1 => ['dump_variation_gene_name_highmem'],
              2 => ['transcript_factory'], 
              3 => ['by_gene_transcript_effect'], 
            },
          },
          { -logic_name => 'dump_variation_gene_name_highmem',
            -module => 'Bio::EnsEMBL::Variation::Pipeline::DumpVariationGeneName',
            -hive_capacity  => 50,
            -parameters => {
              @common_params,
            },
            -rc_name   => 'highmem',
            -flow_into => { 
              2 => ['transcript_factory'], 
              3 => ['by_gene_transcript_effect'], 
            },
          },
          { -logic_name => 'transcript_factory',
            -module => 'Bio::EnsEMBL::Variation::Pipeline::TranscriptFactory',
            -hive_capacity  => 50,
            -analysis_capacity  => 50,
            -parameters => {
              @common_params,
            },
            -rc_name   => 'medmem',
            -flow_into => {
              '2->A' => ['transcript_effect'],
              'A->1' => ['finish_transcript_effect'],
            },
          },
          { -logic_name => 'by_gene_transcript_effect',
            -module => 'Bio::EnsEMBL::Variation::Pipeline::TranscriptEffect',
            -hive_capacity  => 50,
            -analysis_capacity  => 50,
            -parameters => {
              mtmp_table => $self->o('mtmp_table'),
              fasta => $self->o('fasta'),
              disambiguate_single_nucleotide_alleles => $self->o('disambiguate_single_nucleotide_alleles'),
              @common_params,
            },
            -rc_name   => 'default',
          },
          { -logic_name => 'transcript_effect',
            -module => 'Bio::EnsEMBL::Variation::Pipeline::TranscriptEffect',
            -hive_capacity  => 50,
            -analysis_capacity  => 50,
            -parameters => {
              mtmp_table => $self->o('mtmp_table'),
              fasta => $self->o('fasta'),
              disambiguate_single_nucleotide_alleles => $self->o('disambiguate_single_nucleotide_alleles'),
              @common_params,
            },
            -rc_name   => 'medmem',
            -flow_into => {
              -1 => ['transcript_effect_highmem'],
            },
          },
          { -logic_name => 'transcript_effect_highmem',
            -module => 'Bio::EnsEMBL::Variation::Pipeline::TranscriptEffect',
            -hive_capacity  => 50,
            -analysis_capacity  => 50,
            -rc_name => 'highmem',
            -parameters => {
              mtmp_table => $self->o('mtmp_table'),
              fasta => $self->o('fasta'),
              disambiguate_single_nucleotide_alleles => $self->o('disambiguate_single_nucleotide_alleles'),
              @common_params,
            },
          },
          { -logic_name => 'finish_transcript_effect',
            -module => 'Bio::EnsEMBL::Variation::Pipeline::FinishTranscriptEffect',
            -parameters => {
              @common_params,
            },
            -rc_name   => 'default',
          },
          { -logic_name => 'web_index_load',
            -module => 'Bio::EnsEMBL::Variation::Pipeline::LoadWebIndexFiles',
            -parameters => {
              @common_params,
            },
            -rc_name   => 'default',
          },
          { -logic_name => 'rebuild_tv_indexes',
            -module => 'Bio::EnsEMBL::Variation::Pipeline::RebuildIndexes',
            -parameters => {
              @common_params,
            },
            -rc_name   => 'default',
            -wait_for => 'web_index_load',
            -flow_into => {
              1 => ['update_variation_feature'],
            },
          },
          { -logic_name => 'update_variation_feature',
            -module => 'Bio::EnsEMBL::Variation::Pipeline::UpdateVariationFeature',
            -parameters => {
              @common_params,
            },
            -rc_name   => 'default',
            -flow_into => {
              1 => ['check_transcript_variation']
            },
          },
          { -logic_name => 'check_transcript_variation',
            -module => 'Bio::EnsEMBL::Variation::Pipeline::CheckTranscriptVariation',
            -parameters => {
              @common_params,
            },
            -rc_name   => 'default',
          }
        );
    }

   if ($self->o('run_variation_class') && (($self->o('species') !~ /homo_sapiens|human/) || $self->o('human_default_override'))) {

        push @analyses, (

            {   -logic_name     => 'init_variation_class',
                -module         => 'Bio::EnsEMBL::Variation::Pipeline::InitVariationClass',
                -parameters     => {
                    num_chunks  => 50,
                    
                    @common_params,
                },
                -input_ids      => [{}],
                -hive_capacity  => 1,
                -rc_name        => 'default',
                -wait_for       => ( $self->o('run_transcript_effect') ? [ 'update_variation_feature' ] : [] ),
                -flow_into      => {
                    1 => [ 'finish_variation_class' ],
                    2 => [ 'set_variation_class' ],
                },
            },
            
            {   -logic_name     => 'set_variation_class',
                -module         => 'Bio::EnsEMBL::Variation::Pipeline::SetVariationClass',
                -parameters     => {
                    identify_marker_e => $self->o('identify_marker_e'), 
                    @common_params,
                },
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
                -input_ids      => [],
                -hive_capacity  => 1,
                -rc_name        => 'default',
                -wait_for       => [ 'set_variation_class' ],
                -flow_into      => {},
            },

        );
    }

    return \@analyses;
}

1;

