=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2019] EMBL-European Bioinformatics Institute

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
package Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::ReleaseFrequencyDumps_conf;

use strict;
use warnings;

use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
 # All Hive databases configuration files should inherit from HiveGeneric, directly or indirectly
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

    return {
        %{ $self->SUPER::default_options()
        },    # inherit other stuff from the base class

        hive_auto_rebalance_semaphores => 1,
        hive_default_max_retry_count => 0,
        hive_force_init      => 1,
        hive_use_param_stack => 1,
        ensembl_release    => $self->o('ensembl_release'),

        # include or exclude the following species from the dumps, run for a division or all the species on the server
        species => [],
        antispecies => [],
        division    => [],
        run_all     => 1,

        pipeline_name         => $self->o('pipeline_name'),
        homo_sapiens_dump_dir => $self->o('homo_sapiens_dump_dir'),
        vep_cache_dir         => $self->o('vep_cache_dir'),
        assembly              => $self->o('assembly'),
        release               => $self->o('release'),
        step_size             => 500_000,
        overlap               => 500,
                                 # qw/AFR AMR EAS EUR SAS AA EA gnomAD gnomAD_AFR gnomAD_AMR gnomAD_ASJ gnomAD_EAS gnomAD_FIN gnomAD_NFE gnomAD_OTH gnomAD_SAS/;
        af_keys               => ['AFR', 'AMR', 'EAS', 'EUR', 'SAS'],
        output_file_name      => '1000GENOMES-phase_3',
       
        pipeline_wide_analysis_capacity => 25,        

        pipeline_db => {
            -host   => $self->o('hive_db_host'),
            -port   => $self->o('hive_db_port'),
            -user   => $self->o('hive_db_user'),
            -pass   => $self->o('hive_db_password'),            
            -dbname => $ENV{'USER'} . '_' . $self->o('pipeline_name'),
            -driver => 'mysql',
        },
    };
}

sub resource_classes {
    my ($self) = @_;
    return {
        %{$self->SUPER::resource_classes},
        'default' => { 'LSF' => '-q production-rh74 -R"select[mem>1500] rusage[mem=1500]" -M1500'},
        'medium'  => { 'LSF' => '-q production-rh74 -R"select[mem>4500] rusage[mem=4500]" -M4500'},
        'high'    => { 'LSF' => '-q production-rh74 -R"select[mem>8500] rusage[mem=8500]" -M8500'},
    };
}

sub pipeline_analyses {
  my ($self) = @_;
  my @analyses;
  push @analyses, (
      # pre run checks, directories exist etc
      {   -logic_name => 'init_gvf_files',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
          -input_ids     => [{}],
          -parameters => {
            'homo_sapiens_dump_dir' => $self->o('homo_sapiens_dump_dir'),
            'inputcmd'              => 'find #homo_sapiens_dump_dir#/gvf/homo_sapiens -type f -name "homo_sapiens-chr*.gvf.gz" -printf "%f\n"',
          },
          -flow_into  => {
            '2->A' => {'sort_gvf' => {'gvf_file' => '#_0#'}},
            'A->1' => ['init_add_frequencies_gvf'],
          },
      },
      { -logic_name => 'sort_gvf',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
          'homo_sapiens_dump_dir' => $self->o('homo_sapiens_dump_dir'),       
          'cmd'                   => 'zgrep -h -v ^# #homo_sapiens_dump_dir#/gvf/homo_sapiens/#gvf_file# | sort -k4,4n -k5,5n > #homo_sapiens_dump_dir#/gvf/homo_sapiens/sorted_#gvf_file#',
        },
      },
      { -logic_name => 'init_add_frequencies_gvf',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
          -parameters => {
            'homo_sapiens_dump_dir' => $self->o('homo_sapiens_dump_dir'),
            'inputcmd'              => 'find #homo_sapiens_dump_dir#/gvf/homo_sapiens -type f -name "sorted_homo_sapiens-chr*.gvf.gz" -printf "%f\n"',
          },
          -flow_into  => {
            '2->A' => {'add_frequencies_gvf' => {'input_file' => '#_0#'}},
            'A->1' => ['finish_add_frequencies_gvf'],
          },
      },
      { -logic_name => 'add_frequencies_gvf',
        -module => 'Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::DumpFrequencies',
        -parameters => {
          'file_type'             => 'gvf',
          'homo_sapiens_dump_dir' => $self->o('homo_sapiens_dump_dir'),
          'vep_cache_dir'         => $self->o('vep_cache_dir'),
          'release'               => $self->o('release'),
          'assembly'              => $self->o('assembly'),
          'step_size'             => $self->o('step_size'),
          'overlap'               => $self->o('overlap'),
          'af_keys'               => $self->o('af_keys'),
          'output_file_name'      => $self->o('output_file_name'),
        },
      },
      { -logic_name => 'finish_add_frequencies_gvf',
        -module => 'Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::FinishDumpFrequencies',
      }

  );
  return \@analyses;
}
1;
