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
        hive_auto_rebalance_semaphores => 0,  # do not attempt to rebalance semaphores periodically by default
        hive_no_init            => 0, # setting it to 1 will skip pipeline_create_commands (useful for topping up)
        hive_root_dir           => $ENV{'HOME'} . '/bin/ensembl-hive',
        ensembl_cvs_root_dir    => $ENV{'HOME'} . '/bin',
        hive_db_port            => 4521,
        hive_db_user            => 'ensadmin',
        hive_db_host            => 'mysql-ens-var-prod-2.ebi.ac.uk',
        debug                   => 0,
        run_variant_qc          => 1,
        use_fasta_files         => 0,
        flank_seq_length        => 150,
        algn_score_threshold    => 0.95,
        max_map_weight          => 5,
        use_prior_for_filtering => 1,
        map_to_chrom_only       => 1,
        entries_per_file        => 200000,
        mode                    => 'remap_db_table', # options: remap_db_table (default), remap_multi_map, remap_alt_loci, remap_read_coverage, remap_post_projection, remap_svf_post_projection, remap_svf, remap_qtls
        dump_multi_map          => 1,
        feature_table           => 'variation_feature',
        feature_table_failed_projection => 0, #'variation_feature_failed'
        feature_table_projection => 0, #'variation_feature_projection',
        individuals             => '',
        pipeline_dir            => $self->o('pipeline_dir'),
        bam_files               => $self->o('pipeline_dir') . '/bam_files',
        dump_features           => $self->o('pipeline_dir') . '/dump_features',
        fasta_files             => $self->o('pipeline_dir') . '/fasta_files',
        filtered_mappings       => $self->o('pipeline_dir') . '/filtered_mappings',
        load_features           => $self->o('pipeline_dir') . '/load_features',
        mapping_results         => $self->o('pipeline_dir') . '/mapping_results', 
        statistics              => $self->o('pipeline_dir') . '/statistics',   
        dump_mapped_features    => $self->o('pipeline_dir') . '/dump_mapped_features',       
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
        %{$self->SUPER::pipeline_wide_parameters},          # here we inherit anything from the base class
        mode                         => $self->o('mode'),
        registry_file                => $self->o('registry_file'),
        registry_file_newasm         => $self->o('registry_file_newasm'),
        dump_features_dir            => $self->o('dump_features'),
        dump_mapped_features_dir     => $self->o('dump_mapped_features'),
        filtered_mappings_dir        => $self->o('filtered_mappings'),
        load_features_dir            => $self->o('load_features'),
        statistics_dir               => $self->o('statistics'),
        fasta_files_dir              => $self->o('fasta_files'),
        bam_files_dir                => $self->o('bam_files'),
        mapping_results_dir          => $self->o('mapping_results'),
        old_assembly_fasta_file_dir  => $self->o('old_assembly'),
        new_assembly_fasta_file_dir  => $self->o('new_assembly'),
        new_assembly_fasta_file_name => $self->o('new_assembly_file_name'),
        species                      => $self->o('species'),
        bwa_dir                      => $self->o('bwa_location'),
        samtools_dir                 => $self->o('samtools_location'),
        pipeline_dir                 => $self->o('pipeline_dir'),
        flank_seq_length             => $self->o('flank_seq_length'),	
        use_fasta_files              => $self->o('use_fasta_files'),
        feature_table                => $self->o('feature_table'),
        feature_table_projection     => $self->o('feature_table_projection'),
        feature_table_failed_projection => $self->o('feature_table_failed_projection'),
        individuals                  => $self->o('individuals'),
        algn_score_threshold         => $self->o('algn_score_threshold'),
        max_map_weight               => $self->o('max_map_weight'),
        use_prior_for_filtering      => $self->o('use_prior_for_filtering'),
        map_to_chrom_only            => $self->o('map_to_chrom_only'),
        entries_per_file             => $self->o('entries_per_file'),
        debug                        => $self->o('debug'),
        dump_multi_map               => $self->o('dump_multi_map'),
    };
}

sub resource_classes {
    my ($self) = @_;
    return {
        %{$self->SUPER::resource_classes},  # inherit 'default' from the parent class
            'default_mem' => { 'LSF' => '-R"select[mem>2500] rusage[mem=2500]" -M2500'}, 
            'high_mem'    => { 'LSF' => '-R"select[mem>5500] rusage[mem=5500]" -M5500'}, 
    };
}

sub pipeline_analyses {
    my ($self) = @_;
    my @analyses;
    push @analyses, (
        {
            -logic_name => 'pre_run_checks',
            -module     => 'Bio::EnsEMBL::Variation::Pipeline::Remapping::PreRunChecks',
            -input_ids  => [{},],
            -flow_into  => {
                1 => ['init_mapping']
            },
        },
        {   
            -logic_name        => 'init_mapping', 
            -module            => 'Bio::EnsEMBL::Variation::Pipeline::Remapping::InitMapping',
            -rc_name           => 'default_mem',
            -analysis_capacity => 5,
            -flow_into => { 
                '2->A' => ['run_mapping'],
                'A->1' => ['finish_mapping']
            },		
        },
        {
            -logic_name => 'run_mapping',
            -module     => 'Bio::EnsEMBL::Variation::Pipeline::Remapping::RunMapping',
            -rc_name    => 'high_mem',
        },
        {	
            -logic_name => 'finish_mapping',
            -module     => 'Bio::EnsEMBL::Variation::Pipeline::Remapping::FinishMapping',
            -flow_into  => {
                1 => ['init_parse_mapping'],
            },
        },
        {
            -logic_name        => 'init_parse_mapping', 
            -module            => 'Bio::EnsEMBL::Variation::Pipeline::Remapping::InitParseMapping',
            -rc_name           => 'default_mem',
            -analysis_capacity => 5,
            -flow_into => { 
                '2->A' => ['parse_mapping'],
                'A->1' => ['finish_parse_mapping']
            },		
        },
        {
            -logic_name        => 'parse_mapping',
            -module            => 'Bio::EnsEMBL::Variation::Pipeline::Remapping::ParseMapping',
            -analysis_capacity => 5,
            -rc_name           => 'default_mem',
        },
        {
            -logic_name => 'finish_parse_mapping',
            -module     => 'Bio::EnsEMBL::Variation::Pipeline::Remapping::FinishParseMapping',
            -flow_into => {
                1 => ['init_filter_mapping'],
            },
        },
        {
            -logic_name        => 'init_filter_mapping',
            -module            => 'Bio::EnsEMBL::Variation::Pipeline::Remapping::InitFilterMapping',
            -rc_name           => 'default_mem',
            -analysis_capacity => 5,
            -flow_into => {
                '2->A' => ['filter_mapping'],
                'A->1' => ['finish_filter_mapping']
            },
        },
        {
            -logic_name        => 'filter_mapping',
            -module            => 'Bio::EnsEMBL::Variation::Pipeline::Remapping::FilterMapping',
            -analysis_capacity => 5,
            -rc_name           => 'default_mem',
        },
        {
            -logic_name => 'finish_filter_mapping',
            -module     => 'Bio::EnsEMBL::Variation::Pipeline::Remapping::FinishFilterMapping',
            -flow_into => {
                1 => ['load_mapping'],
            },

        },
    );
    if ($self->o('run_variant_qc')) {
      push @analyses, (
        {
            -logic_name => 'load_mapping',
            -module     => 'Bio::EnsEMBL::Variation::Pipeline::Remapping::LoadMapping',
            -flow_into => {
                1 => ['variant_qc'],
            },
        },
        {
            -logic_name => 'variant_qc',
            -module     => 'Bio::EnsEMBL::Variation::Pipeline::Remapping::InitVariantQC',
        },
      );
    } else {
      push @analyses, (
        {
            -logic_name => 'load_mapping',
            -module     => 'Bio::EnsEMBL::Variation::Pipeline::Remapping::LoadMapping',
        },
      );
    }
    return \@analyses;
}

1;

