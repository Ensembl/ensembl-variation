=head1 LICENSE

Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');

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

        # general pipeline options that you should change to suit your environment
       
        hive_force_init => 1,
        hive_use_param_stack => 0,
        hive_use_triggers => 0,
        hive_auto_rebalance_semaphores => 0, 
        hive_no_init => 0,
        # the location of your checkout of the ensembl API (the hive looks for SQL files here)
        
        ensembl_cvs_root_dir    => $ENV{'HOME'} . '/DEV',
        hive_root_dir           => $ENV{'HOME'} . '/DEV/ensembl-hive', 
        # a name for your pipeline (will also be used in the name of the hive database)
        
        pipeline_name           => 'variation_consequence',

        # a directory to keep hive output files and your registry file, you should
        # create this if it doesn't exist

        pipeline_dir            => '/lustre/scratch109/ensembl/at7/release_77/human/' . $self->o('pipeline_name'),

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
        # EBI farm users should either change these here, or override them on the
        # command line to suit the EBI farm. The names of each option hopefully
        # reflect their usage, but you may want to change the details (memory
        # requirements, queue parameters etc.) to suit your own data
        
        default_lsf_options => '-R"select[mem>2000] rusage[mem=2000]" -M2000',
        urgent_lsf_options  => '-q yesterday -R"select[mem>2000] rusage[mem=2000]" -M2000',
        highmem_lsf_options => '-R"select[mem>15000] rusage[mem=15000]" -M15000', # this is Sanger LSF speak for "give me 15GB of memory"
        long_lsf_options    => '-q long -R"select[mem>2000] rusage[mem=2000]" -M2000',

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

        hive_db_host    => 'ens-variation',
        hive_db_port    => 3306,
        hive_db_user    => 'ensadmin',

        pipeline_db => {
            -host   => $self->o('hive_db_host'),
            -port   => $self->o('hive_db_port'),
            -user   => $self->o('hive_db_user'),
            -pass   => $self->o('hive_db_password'),            
            -dbname => $ENV{'USER'}.'_'.$self->o('pipeline_name').'_'.$self->o('species'),
            -driver => 'mysql',
        },
    };
}


sub resource_classes {
    my ($self) = @_;
    return {
          'default' => { 'LSF' => $self->o('default_lsf_options') },
          'urgent'  => { 'LSF' => $self->o('urgent_lsf_options')  },
          'highmem' => { 'LSF' => $self->o('highmem_lsf_options') },
          'long'    => { 'LSF' => $self->o('long_lsf_options')    },
    };
}

sub pipeline_analyses {
    my ($self) = @_;

    my @common_params = (
        ensembl_registry    => $self->o('reg_file'),
        species             => $self->o('species'),
    );
   
    my @analyses;

    if ($self->o('run_transcript_effect')) {

        push @analyses, (
            
            {   -logic_name => 'init_transcript_effect',
                -module     => 'Bio::EnsEMBL::Variation::Pipeline::InitTranscriptEffect',
                -parameters => {
                    include_lrg => $self->o('include_lrg'),
                    limit_biotypes => $self->o('limit_biotypes'),
                    @common_params,
                },
                -input_ids  => [{}],
                -rc_name    => 'default',
                -flow_into  => {
                    2 => [ 'rebuild_tv_indexes' ],
                    3 => [ 'update_variation_feature' ],
                    4 => [ 'transcript_effect' ],
                    5 => [ 'check_transcript_variation' ],

                },
            },

            {   -logic_name     => 'transcript_effect',
                -module         => 'Bio::EnsEMBL::Variation::Pipeline::TranscriptEffect',
                -parameters     => { 
                    disambiguate_single_nucleotide_alleles => $self->o('disambiguate_single_nucleotide_alleles'), 
                    @common_params,
                },
                -input_ids      => [],
                -hive_capacity  => $self->o('transcript_effect_capacity'),
                -rc_name        => 'default',
                -flow_into      => {},
            },

            {   -logic_name     => 'rebuild_tv_indexes',
                -module         => 'Bio::EnsEMBL::Variation::Pipeline::RebuildIndexes',
                -parameters     => {
                    @common_params,
                },
                -input_ids      => [],
                -hive_capacity  => 1,
                -rc_name        => 'urgent',
                -wait_for       => [ 'transcript_effect' ],
                -flow_into      => {},
            },
        
            {   -logic_name     => 'check_transcript_variation',
                -module         => 'Bio::EnsEMBL::Variation::Pipeline::CheckTranscriptVariation',
                -parameters     => {
                    pipeline_dir  => $self->o('pipeline_dir'),
                    @common_params,
                },
                -input_ids      => [],
                -hive_capacity  => 1,
                -rc_name        => 'default',
                -wait_for       => [ 'rebuild_tv_indexes' ],
                -flow_into      => {},
            },

            {   -logic_name     => 'update_variation_feature',
                -module         => 'Bio::EnsEMBL::Variation::Pipeline::UpdateVariationFeature',
                -parameters     => {
                    @common_params,
                },
                -input_ids      => [],
                -hive_capacity  => 1,
                -rc_name        => 'urgent',
                -wait_for       => [ 'rebuild_tv_indexes' ],
                -flow_into      => {},
            }, 

        );
    }

    if ($self->o('run_variation_class')) {

        push @analyses, (

            {   -logic_name     => 'init_variation_class',
                -module         => 'Bio::EnsEMBL::Variation::Pipeline::InitVariationClass',
                -parameters     => {
                    num_chunks  => 50,
                    
                    @common_params,
                },
                -input_ids      => [{}],
                -hive_capacity  => 1,
                -rc_name        => 'highmem',
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
                -rc_name        => 'urgent',
                -wait_for       => [ 'set_variation_class' ],
                -flow_into      => {},
            },

        );
    }

    return \@analyses;
}

1;

