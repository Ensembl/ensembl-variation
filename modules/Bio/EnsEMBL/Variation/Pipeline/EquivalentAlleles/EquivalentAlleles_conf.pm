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

package Bio::EnsEMBL::Variation::Pipeline::EquivalentAlleles::EquivalentAlleles_conf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');

sub default_options {
  my ($self) = @_;

  my $login = `whoami`;
  chomp $login;

  return {

    # general pipeline options that you should change to suit your environment
       
    hive_force_init => 1,
    hive_use_param_stack => 0,
    hive_use_triggers => 0,
    hive_auto_rebalance_semaphores => 0, 
    hive_no_init => 0,

    # the location of your checkout of the ensembl API (the hive looks for SQL files here)
        
    ensembl_cvs_root_dir    => $ENV{'HOME'} . '/bin/',
    hive_root_dir           => $ENV{'HOME'} . '/bin/ensembl-hive', 

    pipeline_name           => 'equivalent_alleles',

    # a directory to keep hive output files and your registry file, you should
    # create this if it doesn't exist

    pipeline_dir            => '/gpfs/nobackup/ensembl/'. $login .'/'. $self->o('pipeline_name') . '/'.  $self->o('species'),

    # a directory where hive workers will dump STDOUT and STDERR for their jobs
        
    output_dir              => $self->o('pipeline_dir').'/hive_output',

    # a standard ensembl registry file containing connection parameters
    # for your target database(s) 
      
    reg_file                => $self->o('pipeline_dir').'/ensembl.registry',

    # configuration for the various resource options used in the pipeline
    # EBI farm users should either change these here, or override them on the
    # command line to suit the EBI farm.
        
    default_lsf_options => '-R"select[mem>2000] rusage[mem=2000]" -M2000',
    medium_lsf_options  => '-R"select[mem>4000] rusage[mem=4000]" -M4000',


    # size of region to be checked in a single job
    region_size =>  100000000,

    # size of bin to be checked for equivalent alleles
    bin_size  =>  5000000,

    ## overlap between checking bins
    overlap   => 1000,

    # number of workers used for the parallelisable analysis
    capacity  => 30,


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
        -dbname => $ENV{'USER'}.'_'.$self->o('pipeline_name') . '_' . $self->o('species'),
        -driver => 'mysql',
    },
  };
}

sub resource_classes {
  my ($self) = @_;
  return {
      'default' => { 'LSF' => $self->o('default_lsf_options') },
      'medium'  => { 'LSF' => $self->o('medium_lsf_options')  },
  };
}

sub pipeline_analyses {
  my ($self) = @_;

  my @common_params = (
      ensembl_registry    => $self->o('reg_file'),
      species             => $self->o('species'),
      pipeline_dir        => $self->o('pipeline_dir'),
  );
   
   
  my @analyses;
 
  push @analyses, (
    
    { -logic_name => 'init_equivalent_alleles',
      -module     => 'Bio::EnsEMBL::Variation::Pipeline::EquivalentAlleles::InitEquivalentAlleles',
      -parameters => {
          region_size       => $self->o('region_size'),
          overlap           => $self->o('overlap'),
          @common_params,
      },
      -input_ids  => [{}],
      -hive_capacity  => -1,
      -rc_name    => 'medium',
      -flow_into  => {
           2 => [ 'find_equivalent_alleles' ],
           3 => [ 'finish_equivalent_alleles' ]
      }
    },
    {  -logic_name     => 'find_equivalent_alleles',
       -module         => 'Bio::EnsEMBL::Variation::Pipeline::EquivalentAlleles::FindEquivalentAlleles',        
       -parameters     => {   
             overlap  => $self->o('overlap'),
             bin_size => $self->o('bin_size'),
             @common_params,
        },
       -input_ids        => [],
       -hive_capacity    => $self->o('capacity'),
       -max_retry_count  => 0,
       -rc_name          => 'default',
       -wait_for         => [ 'init_equivalent_alleles' ],
       -flow_into        => {},
    },

    {  -logic_name     => 'finish_equivalent_alleles',
       -module         => 'Bio::EnsEMBL::Variation::Pipeline::EquivalentAlleles::FinishEquivalentAlleles', 
       -parameters     => {
           @common_params,
        },
        -input_ids      => [],
        -hive_capacity  => -1,
        -rc_name        => 'default',
        -wait_for       => [ 'init_equivalent_alleles','find_equivalent_alleles' ],
        -flow_into      => {},
    },                   
   
  );
 
    return \@analyses;
}

1;

