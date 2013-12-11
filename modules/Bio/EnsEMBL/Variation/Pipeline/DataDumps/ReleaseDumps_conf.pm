=head1 LICENSE

 Copyright (c) 1999-2013 The European Bioinformatics Institute and
 Genome Research Limited.  All rights reserved.

 This software is distributed under a modified Apache license.
 For license details, please see

   http://www.ensembl.org/info/about/legal/code_licence.html

=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <dev@ensembl.org>.

 Questions may also be sent to the Ensembl help desk at
 <helpdesk@ensembl.org>.

=cut

package Bio::EnsEMBL::Variation::Pipeline::DataDumps::ReleaseDumps_conf;

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
		hive_force_init      => 1,
		hive_use_param_stack => 0,
        hive_use_triggers    => 0,
        # general pipeline options that you should change to suit your environment
		hive_root_dir        => $ENV{'HOME'} . '/HEAD/ensembl-hive',	
        ensembl_cvs_root_dir => $ENV{'HOME'} . '/HEAD',


        # a name for your pipeline (will also be used in the name of the hive database)
        pipeline_name           => 'test_data_dumps_e' . $self->o('ensembl_release'),

        # a directory to keep hive output files and your registry file, you should
        # create this if it doesn't exist

        pipeline_dir            => '/lustre/scratch110/ensembl/at7/develop/test/data_dumps/',

        # a directory where hive workers will dump STDOUT and STDERR for their jobs
        # if you use lots of workers this directory can get quite big, so it's
        # a good idea to keep it on lustre, or some other place where you have a 
        # healthy quota!
        output_dir              => $self->o('pipeline_dir') . '/hive_output',

        # a standard ensembl registry file containing connection parameters
        # for your target database(s) (and also possibly aliases for your species
        # of interest that you can then supply to init_pipeline.pl with the -species
        # option)
        
        reg_file      => $self->o('pipeline_dir')         . '/ensembl.registry.75',
		config_file   => $self->o('pipeline_dir')         . '/ovis_aries_75.json',
		tmp_dir       => $self->o('pipeline_dir')         . '/tmp', 
        script_dir    => $self->o('ensembl_cvs_root_dir') . '/ensembl-variation/scripts',
		gvf_validator => '/nfs/users/nfs_a/at7/tools/gvf_validator',
		so_file       => '/nfs/users/nfs_a/at7/obo/e74/so.obo',
		debug         => 0, 

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
            -dbname => $ENV{'USER'} . '_' . $self->o('pipeline_name'),
			-driver => 'mysql',
        },
    };
}

sub pipeline_wide_parameters {
    my ($self) = @_;
    return {
        %{$self->SUPER::pipeline_wide_parameters}, # here we inherit anything from the base class
		release          => $self->o('ensembl_release'),
        ensembl_registry => $self->o('reg_file'),
        config_file      => $self->o('config_file'),
        output_path      => $self->o('pipeline_dir'),
        debug            => $self->o('debug'),
        data_dump_dir    => $self->o('pipeline_dir'),
		tmp_dir          => $self->o('tmp_dir'),		
        script_dir       => $self->o('script_dir'), 
		so_file          => $self->o('so_file'),	
		gvf_validator    => $self->o('gvf_validator'),
		config_file      => $self->o('config_file'),
    };
}

sub resource_classes {
	my ($self) = @_;
	return {
		%{$self->SUPER::resource_classes},
	    'default' => { 'LSF' => '-R"select[mem>4000] rusage[mem=4000]" -M4000'},
        'urgent'  => { 'LSF' => '-q yesterday -R"select[mem>2000] rusage[mem=2000]" -M2000'},
        'highmem' => { 'LSF' => '-R"select[mem>15000] rusage[mem=15000]" -M15000'}, # this is Sanger LSF speak for "give me 15GB of memory"
        'long'    => { 'LSF' => '-q long -R"select[mem>2000] rusage[mem=2000]" -M2000'},
	};
}

sub pipeline_analyses {
    my ($self) = @_;
    my @analyses;
		push @analyses, (
			{	-logic_name => 'pre_run_checks_gvf_dumps',
				-module => 'Bio::EnsEMBL::Variation::Pipeline::DataDumps::PreRunChecks',
				-input_ids => [{},],
				-max_retry_count => 1,
				-flow_into => {
					1 => ['species_factory_gvf_dumps'],
				},
                -parameters => {
					'file_type' => 'gvf',	
				},
			},
			{	-logic_name => 'species_factory_gvf_dumps',
				-module => 'Bio::EnsEMBL::Variation::Pipeline::DataDumps::SpeciesFactory', 
				-analysis_capacity => 3,
				-flow_into => { 
					'2->A' => ['init_dump'],
					'A->1' => ['report_gvf_dumps']			
				 },		
			},
            {   -logic_name => 'init_dump',
                -module => 'Bio::EnsEMBL::Variation::Pipeline::DataDumps::InitSubmitJob',
				-analysis_capacity => 5,
				-flow_into => {
					1 => ['submit_job_gvf_dumps'],
				},
                -parameters => {
					'file_type' => 'gvf',
					'job_type'  => 'dump',
				},
				-rc_name => 'default',
            },
			{	-logic_name => 'submit_job_gvf_dumps',
				-module => 'Bio::EnsEMBL::Variation::Pipeline::DataDumps::SubmitJob',
				-max_retry_count => 1,
				-rc_name => 'default',
			},
			{	-logic_name => 'report_gvf_dumps',
				-module => 'Bio::EnsEMBL::Variation::Pipeline::DataDumps::Report',
                -parameters => {
					'file_type'   => 'gvf',
					'report_name' => 'ReportDumpGVF',
				},
				-flow_into => {1 => ['species_factory_validate_gvf']},
			},
			{	-logic_name => 'species_factory_validate_gvf',	
				-module => 'Bio::EnsEMBL::Variation::Pipeline::DataDumps::SpeciesFactory', 
				-analysis_capacity => 3,
				-flow_into => { 
					'2->A' => ['join_dump'],
					'A->1' => ['summary_validate_gvf']			
				 },		
                -parameters => {
				    'file_type' => 'gvf',		
				},
			},
			{	-logic_name => 'join_dump',
				-module => 'Bio::EnsEMBL::Variation::Pipeline::DataDumps::JoinDump',
				-flow_into => {
					1 => ['validate_gvf'],
				}
			},
			{	-logic_name => 'validate_gvf',
			 	-module => 'Bio::EnsEMBL::Variation::Pipeline::DataDumps::ValidateDump',
                -parameters => {
				    'file_type'     => 'gvf',		
					'so_file'       => $self->o('so_file'),	
					'gvf_validator' => $self->o('gvf_validator'),
				},
			},
			{	-logic_name => 'summary_validate_gvf', # die if Error?
			 	-module => 'Bio::EnsEMBL::Variation::Pipeline::DataDumps::SummariseValidation',
                -parameters => {
				    'file_type' => 'gvf',		
				},
				-flow_into => {
					1 => ['clean_up_gvf_dumps']
				},	
			},
			{	-logic_name => 'clean_up_gvf_dumps',
			 	-module => 'Bio::EnsEMBL::Variation::Pipeline::DataDumps::CleanUp',
				-flow_into => {
					1 => ['pre_run_checks_gvf2vcf']
				},	
				-parameters => {
					'file_type' => 'gvf',
				},
			},
			{	-logic_name => 'pre_run_checks_gvf2vcf', 
				-module => 'Bio::EnsEMBL::Variation::Pipeline::DataDumps::PreRunChecks',
				-max_retry_count => 1,
				-flow_into => {
					1 => ['species_factory_gvf2vcf'],
				},
                -parameters => {
					'file_type' => 'vcf',
				},
			},
			{	-logic_name => 'species_factory_gvf2vcf',
				-module => 'Bio::EnsEMBL::Variation::Pipeline::DataDumps::SpeciesFactory', 
				-analysis_capacity => 3,
				-flow_into => { 
					'2->A' => ['init_parse'],
					'A->1' => ['summary_validate_vcf']			
				 },		
                -parameters => { 
					'file_type' => 'vcf',		
				},
			},
			{	-logic_name => 'init_parse',
				-module => 'Bio::EnsEMBL::Variation::Pipeline::DataDumps::InitSubmitJob',
				-flow_into => {
					1 => ['submit_job_gvf2vcf'],
				},
                -parameters => {
					'file_type' => 'vcf',
					'job_type'  => 'parse',
				},
			},
			{	-logic_name => 'submit_job_gvf2vcf',
				-module => 'Bio::EnsEMBL::Variation::Pipeline::DataDumps::SubmitJob',
				-flow_into => {
					1 => ['validate_vcf'],	
				},
			},
			{	-logic_name => 'validate_vcf',
				-module => 'Bio::EnsEMBL::Variation::Pipeline::DataDumps::ValidateVCF',
				-parameters => {
					'file_type' => 'vcf',
				},
			}, 
			{	-logic_name => 'summary_validate_vcf', # die if Error?
			 	-module => 'Bio::EnsEMBL::Variation::Pipeline::DataDumps::SummariseValidation',
				-parameters => {
					'file_type' => 'vcf',
				},
				-flow_into => {
					1 => ['clean_up_vcf_dumps'],
				},
			},
			{	-logic_name => 'clean_up_vcf_dumps',
			 	-module => 'Bio::EnsEMBL::Variation::Pipeline::DataDumps::CleanUp',
				-parameters => {
					'file_type' => 'vcf',
				},
				-flow_into => {
					1 => ['finish'],
				},
			},
			{	-logic_name => 'finish', # readme, 
				-module     => 'Bio::EnsEMBL::Variation::Pipeline::DataDumps::FinishDump',
			},
		);
    return \@analyses;
}

1;

