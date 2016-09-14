=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

package Bio::EnsEMBL::Variation::Pipeline::TranscriptHaplotypes::TranscriptHaplotypes_conf;

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
    # this pipeline requires you have ensembl-variation and ensembl-tools in this dir as
    ensembl_cvs_root_dir    => $ENV{'HOME'} . '/Git',
    hive_root_dir           => $ENV{'HOME'} . '/Git/ensembl-hive', 
    
    # a name for your pipeline (will also be used in the name of the hive database)    
    pipeline_name           => 'transcript_haplotypes_'.$self->o('species'),

    # a directory to keep hive output files and your registry file, you should
    # create this if it doesn't exist
    pipeline_dir            => '/lustre/scratch109/ensembl/'.$ENV{'USER'}.'/'.$self->o('pipeline_name').'/'.$self->o('species'),

    # a directory where hive workers will dump STDOUT and STDERR for their jobs
    # if you use lots of workers this directory can get quite big, so it's
    # a good idea to keep it on lustre, or some other place where you have a 
    # healthy quota!
    output_dir              => $self->o('pipeline_dir').'/hive_output',

    # a standard ensembl registry file containing connection parameters
    # for your target database(s) (and also possibly aliases for your species
    # of interest that you can then supply to init_pipeline.pl with the -species
    # option)
    
    ensembl_registry        => $self->o('pipeline_dir').'/ensembl.registry',

    # path to a JSON config file for genotype VCFs
    json_config             => undef,

    # configuration for the various resource options used in the pipeline
    # EBI farm users should either change these here, or override them on the
    # command line to suit the EBI farm. The names of each option hopefully
    # reflect their usage, but you may want to change the details (memory
    # requirements, queue parameters etc.) to suit your own data
        
    default_lsf_options => '-R"select[mem>4000] rusage[mem=4000]" -M4000',
    urgent_lsf_options  => '-q yesterday -R"select[mem>2000] rusage[mem=2000]" -M2000',
    highmem_lsf_options => '-q basement -R"select[mem>15000] rusage[mem=15000]" -M15000', # this is Sanger LSF speak for "give me 15GB of memory"
    long_lsf_options    => '-q long -R"select[mem>4000] rusage[mem=4000]" -M4000',

    # init_pipeline.pl will create the hive database on this machine, naming it
    # <username>_<pipeline_name>, and will drop any existing database with this
    # name

    hive_db_host    => 'ens-variation3',
    hive_db_port    => 3306,
    hive_db_user    => 'ensadmin',

    pipeline_db => {
      -host   => $self->o('hive_db_host'),
      -port   => $self->o('hive_db_port'),
      -user   => $self->o('hive_db_user'),
      -pass   => $self->o('hive_db_password'),            
      -dbname => $ENV{'USER'}.'_'.$self->o('pipeline_name'),
      -driver => 'mysql',
    },
    
    debug => 0,
    qc => 1,
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

  my @common_params = map {$_ => $self->o($_) || undef} qw(
    ensembl_registry
    species
    pipeline_dir
    json_config
  );
   
  my @analyses = (
    {
      -logic_name    => 'init_transcript_haplotypes',
      -module        => 'Bio::EnsEMBL::Variation::Pipeline::TranscriptHaplotypes::InitTranscriptHaplotypes',
      -parameters    => {
        @common_params
      },
      -input_ids     => [{}],
      -rc_name       => 'long',
      -hive_capacity => 1,
      -flow_into     => {
        '1' => ['dump_transcript_haplotypes'],
        '2' => ['dump_transcript_haplotypes_highmem'],
        '3' => ['finish_transcript_haplotypes']
      },
    },
    {
      -logic_name    => 'dump_transcript_haplotypes',
      -module        => 'Bio::EnsEMBL::Variation::Pipeline::TranscriptHaplotypes::DumpTranscriptHaplotypes',
      -parameters    => {
        @common_params,
        'filter_frequency'
      },
      -rc_name       => 'default',
      -analysis_capacity => 30,
      -flow_into      => {
        -1 => ['dump_transcript_haplotypes_highmem'],
      }
    },
    {
      -logic_name    => 'dump_transcript_haplotypes_highmem',
      -module        => 'Bio::EnsEMBL::Variation::Pipeline::TranscriptHaplotypes::DumpTranscriptHaplotypes',
      -parameters    => {
        @common_params,
        'filter_frequency'
      },
      -rc_name       => 'highmem',
      -analysis_capacity => 20,
      -can_be_empty   => 1,
    },
    {
      -logic_name    => 'finish_transcript_haplotypes',
      -module        => 'Bio::EnsEMBL::Variation::Pipeline::TranscriptHaplotypes::FinishTranscriptHaplotypes',
      -parameters    => {
        @common_params
      },
      -rc_name       => 'default',
      -analysis_capacity => 1,
      -wait_for      => ['dump_transcript_haplotypes', 'dump_transcript_haplotypes_highmem'],
    },
  );

  return \@analyses;
}

1;

