package Bio::EnsEMBL::Variation::Pipeline::BackupDB::Backup_conf;

use strict;
use warnings;
use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');

sub default_options {
    my ($self) = @_;

# The hash returned from this function is used to configure the
# pipeline, you can supply any of these options on the command
# line to override these default values.

#################### MANDATORY OPTIONS ###########################
# -species:          e.g. human                                  #
# -type:             type of DB dump (e.g. 3)                    #
# -hive_db_password: Hive DB password                            #
# -db_host:          DB host (e.g. ensembldb.ensembl.org)        #
# -db_name:          DB name (e.g. homo_sapiens_variation_79_38) #
# -db_user:          DB user (e.g. anonymous)                    #
# -db_port:          DB port (e.g. 3306)                         #
# -db_version:       Ensembl release (e.g. 79)                   #
##################################################################
# e.g.:
# init_pipeline.pl Bio::EnsEMBL::Variation::Pipeline::BackupDB::Backup_conf -species human -type 3 -db_host ensembldb.ensembl.org -db_name homo_sapiens_variation_79_38 -db_user anonymous -db_port 3306 -db_version 79 -hive_db_password XXXXX 

  return {
  
        hive_auto_rebalance_semaphores => 0,
    
        hive_force_init        => 1,
        hive_no_init           => 0,
        hive_use_param_stack   => 0,
        hive_use_triggers      => 0,
        hive_root_dir          => $ENV{'HOME'} . '/DEV/ensembl-hive', # To be updated!
        hive_db_host           => 'ens-variation',
        hive_db_port           => 3306,
        hive_db_user           => 'ensadmin',
        debug                  => 0,

        pipeline_name          => 'automated_db_backup',
        
        # a directory to keep hive output files, you should create this if it doesn't exist
        pipeline_dir           => '/lustre/scratch110/ensembl/'.$ENV{'USER'}.'/'.$self->o('pipeline_name').'/'.$self->o('species'),

        db_driver              => 'mysql',

        src_db_conn            => $self->o('db_driver').'://'.$self->o('db_user').'@'.$self->o('db_host').':'.$self->o('db_port').'/'.$self->o('db_name'),
        
        exclude_tables         => 'MTMP_%',

        date                   => $self->date(),
        
        backup_dir             => $self->o('pipeline_dir').'/'.$self->o('db_version'), # Should be updated!
        output_file_prefix     => $self->o('date').'_'.$self->o('species').'_',
        output_file_suffix     => '_dump_',
        output_extension       => 'sql.gz',
        info_updates           => $self->o('species').'_info_updates.txt',

        output_dir             => $self->o('pipeline_dir').'/hive_output',
        
        # these flags control which parts of the pipeline are run

        run_extract_xml_files  => 1,
        
        small_lsf_options   => '-R"select[mem>1500] rusage[mem=1500]" -M1500',
        default_lsf_options => '-R"select[mem>2000] rusage[mem=2000]" -M2000',
        highmem_lsf_options => '-R"select[mem>15000] rusage[mem=15000]" -M15000', # this is Sanger LSF speak for "give me 15GB of memory"

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

sub backup_type {
  my ($self,$type) = @_;
  my %types_list = (
           1 => 'all',
           2 => 'update',
           3 => 'variation variation_feature',
           4 => 'variarion variation_feature transcript_variation',
     );
  return $types_list{$type};
}

sub resource_classes {
    my ($self) = @_;
    return {
          'small'   => { 'LSF' => $self->o('small_lsf_options')   },
          'default' => { 'LSF' => $self->o('default_lsf_options') },
          'highmem' => { 'LSF' => $self->o('highmem_lsf_options') },
    };
}

sub date {
    my $self = shift;
    my @time = localtime(time());
    
    $time[4]++;

    # add leading zeroes as required
    for my $i(0..4) {
      $time[$i] = "0".$time[$i] if $time[$i] < 10;
    }

    my $time = ($time[5] + 1900)."-".$time[4]."-".$time[3];

    return $time;
}

sub pipeline_analyses {
    my ($self) = @_;
    my @analyses;
    
    push @analyses, (
        {   
            -logic_name    => 'init_backup', 
            -module        => 'Bio::EnsEMBL::Variation::Pipeline::BackupDB::InitBackup',
            -rc_name       => 'small',
            -parameters    => {
               species            => $self->o('species'),
               backup_dir         => $self->o('backup_dir'),
               src_db_conn        => $self->o('src_db_conn'),
               output_file_prefix => $self->o('output_file_prefix'),
               output_file_suffix => $self->o('output_file_suffix'),
               output_extension   => $self->o('output_extension'),
               tables             => $self->backup_type($self->o('type')),
               exclude_tables     => $self->o('exclude_tables'),
               db_name            => $self->o('db_name'),
               info_updates       => $self->o('info_updates'),
            },
            -input_ids     => [{}],
            -wait_for      => [],
            -flow_into     => { 
               '2->A' => ['backup_database'],
               'A->1' => ['finish_backup']
            },		
        },
        {   
            -logic_name        => 'backup_database', 
            -module            => 'Bio::EnsEMBL::Hive::RunnableDB::DatabaseDumper',
            -rc_name           => 'default',
            -input_ids         => [],
            -hive_capacity     => 4,
            -analysis_capacity => 4,
            -wait_for          => [ 'init_backup' ],
            -flow_into         => {},
        },
        {
            -logic_name    => 'finish_backup', 
            -module        => 'Bio::EnsEMBL::Variation::Pipeline::BackupDB::FinishBackup',
            -rc_name       => 'small',
            -parameters    => {
               species            => $self->o('species'),
               backup_dir         => $self->o('backup_dir'),
               src_db_conn        => $self->o('src_db_conn'),
               db_name            => $self->o('db_name'),
               output_file_prefix => $self->o('output_file_prefix'),
               output_file_suffix => $self->o('output_file_suffix'),
               output_extension   => $self->o('output_extension'),
               exclude_tables     => $self->o('exclude_tables'),
               info_updates       => $self->o('info_updates'),
            },
            -input_ids     => [],
            -wait_for      => [ 'backup_database' ],
            -flow_into     => {},
        },
       
    );
    return \@analyses;
}

1;

