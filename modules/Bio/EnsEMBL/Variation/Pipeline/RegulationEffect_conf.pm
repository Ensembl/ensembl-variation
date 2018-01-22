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

package Bio::EnsEMBL::Variation::Pipeline::RegulationEffect_conf;

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
        
        pipeline_name           => 'regulation_effect',

        species                 => $self->o('species'),
        # a directory to keep hive output files and your registry file, you should
        # create this if it doesn't exist
        
        pipeline_dir            => $self->o('pipeline_dir'),

        # a directory where hive workers will dump STDOUT and STDERR for their jobs
        # if you use lots of workers this directory can get quite big, so it's
        # a good idea to keep it on lustre, or some other place where you have a 
        # healthy quota!
        
        output_dir              => $self->o('pipeline_dir') . '/hive_output',

        # a standard ensembl registry file containing connection parameters
        # for your target database(s) (and also possibly aliases for your species
        # of interest that you can then supply to init_pipeline.pl with the -species
        # option)
        
        registry_file           => $self->o('pipeline_dir') . '/ensembl.registry',

        # if set to 1 this option tells the transcript_effect analysis to disambiguate
        # ambiguity codes in single nucleotide alleles, so e.g. an allele string like
        # 'T/M' will be treated as if it were 'T/A/C' (this was a request from ensembl
        # genomes and we don't use it by default in the ensembl variation pipeline)
        
        disambiguate_single_nucleotide_alleles => 0,


        # set this flag to 1 to include LRG transcripts in the transcript effect analysis

        include_lrg => 1,
        include_external_features => 0,
        debug => 0,

        only_update_vf => 0,
        update_vf => 1,
        only_motif_feature => 0,
        only_regulatory_feature => 0,

        # connection parameters for the hive database, you should supply the hive_db_password
        # option on the command line to init_pipeline.pl (parameters for the target database
        # should be set in the registry file defined above)


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
            -dbname => $ENV{'USER'} . '_' . $self->o('pipeline_name') . '_' . $self->o('species'),
            -driver => 'mysql',
        },
    };
}


sub resource_classes {
    my ($self) = @_;
    # configuration for the various resource options used in the pipeline
    # EBI farm users should either change these here, or override them on the
    # command line to suit the EBI farm. The names of each option hopefully
    # reflect their usage, but you may want to change the details (memory
    # requirements, queue parameters etc.) to suit your own data

    return {
          'default' => { 'LSF' => '-qproduction-rh7 -R"select[mem>2000] rusage[mem=2000]" -M2000'},
          'urgent'  => { 'LSF' => '-qproduction-rh7 -R"select[mem>2000] rusage[mem=2000]" -M2000'},
          'highmem' => { 'LSF' => '-qproduction-rh7 -R"select[mem>15000] rusage[mem=15000]" -M15000'},
          'long'    => { 'LSF' => '-qproduction-rh7 -R"select[mem>2000] rusage[mem=2000]" -M2000'},
    };
}

sub pipeline_wide_parameters {
    my ($self) = @_;
    return {
        %{$self->SUPER::pipeline_wide_parameters},
        ensembl_registry => $self->o('registry_file'),
        include_external_features => $self->o('include_external_features'),
        disambiguate_single_nucleotide_alleles => $self->o('disambiguate_single_nucleotide_alleles'),
        debug => $self->o('debug'),
        only_update_vf => $self->o('only_update_vf'),
        update_vf => $self->o('update_vf'),
        only_motif_feature => $self->o('only_motif_feature'),   
        only_regulatory_feature => $self->o('only_regulatory_feature'),
        species => $self->o('species'),
    };
}

sub pipeline_analyses {
    my ($self) = @_;
   
    my @analyses;
    if ($self->o('only_update_vf')) {
        push @analyses, (
            {   -logic_name => 'finish_regulation_effect',
                -module => 'Bio::EnsEMBL::Variation::Pipeline::FinishRegulationEffect',
                -input_ids => [{},],
                -hive_capacity => 1,
            },
        );
    } else {
        push @analyses, (
            {   -logic_name => 'init_regulation_effect',
                -module => 'Bio::EnsEMBL::Variation::Pipeline::InitRegulationEffect',
                -hive_capacity => 1,
                -rc_name => 'highmem',
                -input_ids => [{},],
                -flow_into => {
                    '2->A' => ['regulation_effect'],
                    'A->1' => ['finish_regulation_effect'],
                },
            },
            {   -logic_name => 'regulation_effect',
                -module => 'Bio::EnsEMBL::Variation::Pipeline::RegulationEffect',
                -rc_name => 'default',
                -hive_capacity  =>  20,
            }, 
            {   -logic_name => 'finish_regulation_effect',
                -module => 'Bio::EnsEMBL::Variation::Pipeline::FinishRegulationEffect',
                -hive_capacity => 1,
            },
        );
    }
    return \@analyses;
}

1;
