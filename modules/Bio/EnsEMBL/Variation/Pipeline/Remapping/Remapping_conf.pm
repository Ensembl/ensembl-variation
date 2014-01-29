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

package Bio::EnsEMBL::Variation::Pipeline::Remapping::Remapping_conf;

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
		hive_force_init         => 1,
		hive_use_param_stack    => 0,
		hive_use_triggers       => 0,
		hive_root_dir           => $ENV{'HOME'} . '/HEAD/ensembl-hive',
        ensembl_cvs_root_dir    => $ENV{'HOME'} . '/HEAD',

		debug                   => 1,		
		species                 => 'rattus_norvegicus',
        pipeline_name           => 'remapping',

        pipeline_dir            => '/lustre/scratch109/ensembl/at7/remapping/rattus_norvegicus/',

        output_dir              => $self->o('pipeline_dir') . '/hive_output',

        reg_file                     => $self->o('pipeline_dir' ). '/ensembl.registry',
		fasta_files_dir              => $self->o('pipeline_dir') . '/fasta_files/',	
		bam_files_dir                => $self->o('pipeline_dir') . '/bam_files/',
		old_assembly_fasta_file_dir  => $self->o('pipeline_dir') . '/old_assembly/',
		new_assembly_fasta_file_dir  => $self->o('pipeline_dir') . '/new_assembly/',
		new_assembly_fasta_file_name => 'Rattus_norvegicus.Rnor_5.0.74.dna.toplevel.fa',
	    mapping_results_dir          => $self->o('pipeline_dir') . '/mapping_results/',

		tool_dir => '/nfs/users/nfs_a/at7/tools/',
		bwa_dir => 'bwa-0.7.5a',
		samtools_dir => 'samtools-0.1.19',

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
            -dbname => $ENV{'USER'} . '_' . $self->o('pipeline_name') . '_' . $self->o('species'),
			-driver => 'mysql',
        },
    };
}

sub pipeline_wide_parameters {
    my ($self) = @_;
    return {
        %{$self->SUPER::pipeline_wide_parameters},          # here we inherit anything from the base class
        registry_file                => $self->o('reg_file'),
		fasta_files_dir              => $self->o('fasta_files_dir'),
        bam_files_dir                => $self->o('bam_files_dir'),
		mapping_results_dir          => $self->o('mapping_results_dir'),
		old_assembly_fasta_file_dir  => $self->o('old_assembly_fasta_file_dir'),
		new_assembly_fasta_file_dir  => $self->o('new_assembly_fasta_file_dir'),
		new_assembly_fasta_file_name => $self->o('new_assembly_fasta_file_name'),
		species                      => $self->o('species'),
		tool_dir                     => $self->o('tool_dir'),
		bwa_dir                      => $self->o('bwa_dir'),
		samtools_dir                 => $self->o('samtools_dir'),
		pipeline_dir                 => $self->o('pipeline_dir'),
		flank_seq_length             => 200,	
		generate_fasta_files         => 1,
		debug                        => $self->o('debug'),
    };
}

sub resource_classes {
    my ($self) = @_;
    return {
    	%{$self->SUPER::resource_classes},  # inherit 'default' from the parent class
		'init_mapping'  => { 'LSF' => '-R"select[mem>2500] rusage[mem=2500]" -M2500'}, 
        'run_mapping'   => { 'LSF' => '-R"select[mem>5500] rusage[mem=5500]" -M5500'}, 
		'parse_mapping' => { 'LSF' => '-R"select[mem>2500] rusage[mem=2500]" -M2500'}, 
   };
}

sub pipeline_analyses {
    my ($self) = @_;
   
    my @analyses;
    push @analyses, (
			{	-logic_name => 'pre_run_checks',
				-module => 'Bio::EnsEMBL::Variation::Pipeline::Remapping::PreRunChecks',
				-input_ids => [{},],
				-flow_into => {
					1 => ['init_mapping']
				},
			},
            {   -logic_name => 'init_mapping', 
                -module => 'Bio::EnsEMBL::Variation::Pipeline::Remapping::InitMapping',
				-rc_name => 'init_mapping',
				-analysis_capacity => 5,
                -flow_into => { 
					'2->A' => ['run_mapping'],
					'A->1' => ['finish_mapping']
				},		
            },
			{	-logic_name => 'run_mapping',
				-module => 'Bio::EnsEMBL::Variation::Pipeline::Remapping::RunMapping',
				-rc_name => 'run_mapping',
			},
			{	-logic_name => 'finish_mapping',
				-module => 'Bio::EnsEMBL::Variation::Pipeline::Remapping::FinishMapping',
				-flow_into => {
					1 => ['init_parse_mapping'],
				},
			},
            {   -logic_name => 'init_parse_mapping', 
                -module => 'Bio::EnsEMBL::Variation::Pipeline::Remapping::InitParseMapping',
				-rc_name => 'init_mapping',
				-analysis_capacity => 5,
                -flow_into => { 
					'2->A' => ['parse_mapping'],
					'A->1' => ['finish_parse_mapping']
				},		
            },
            {   -logic_name => 'parse_mapping',
                -module => 'Bio::EnsEMBL::Variation::Pipeline::Remapping::ParseMapping',
				-analysis_capacity => 5,
				-rc_name => 'parse_mapping',
            },
			{	-logic_name => 'finish_parse_mapping',
				-module => 'Bio::EnsEMBL::Variation::Pipeline::Remapping::FinishParseMapping',
			},
        );
    return \@analyses;
}

1;

