=head1 LICENSE
Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2021] EMBL-European Bioinformatics Institute
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

package Bio::EnsEMBL::Variation::Pipeline::Remapping::RemappingStructuralVariationFeature_conf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Variation::Pipeline::Remapping::Remapping_conf');

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
      %{ $self->SUPER::default_options() },
      mode                    => 'remap_svf',
      debug                   => 0,
      debug_sequence_name     => 3,
      flank_seq_length        => 150,
      algn_score_threshold    => 0.95,
      max_map_weight          => 5,
      use_prior_for_filtering => 1,
      map_to_chrom_only       => 1,
      entries_per_file        => 200000,
      feature_table           => 'structural_variation_feature',
      feature_table_mapping_results => 'structural_variation_feature_mapping_results',
      qc_mapped_features      => $self->o('pipeline_dir') . '/qc_mapped_features',
      hive_db_name            => $ENV{'USER'} . '_ehive_remapping_sv_' . $self->o('ensembl_release') . '_' . $self->o('assembly') . '_' . $self->o('species'),
    };
}

sub pipeline_wide_parameters {
    my ($self) = @_;
    return {
      %{$self->SUPER::pipeline_wide_parameters},          # here we inherit anything from the base class
      mode                          => $self->o('mode'),
      flank_seq_length              => $self->o('flank_seq_length'), 
      feature_table                 => $self->o('feature_table'),
      feature_table_mapping_results => $self->o('feature_table_mapping_results'),
      algn_score_threshold          => $self->o('algn_score_threshold'),
      max_map_weight                => $self->o('max_map_weight'),
      use_prior_for_filtering       => $self->o('use_prior_for_filtering'),
      map_to_chrom_only             => $self->o('map_to_chrom_only'),
      entries_per_file              => $self->o('entries_per_file'),
      debug                         => $self->o('debug'),
      debug_sequence_name           => $self->o('debug_sequence_name'),
      qc_mapped_features_dir        => $self->o('qc_mapped_features'),
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
            -max_retry_count  => 0,
            -input_ids  => [{},],
            -flow_into  => {
                1 => ['init_mapping']
            },
        },
        {   
            -logic_name        => 'init_mapping', 
            -module            => 'Bio::EnsEMBL::Variation::Pipeline::Remapping::InitStructuralVariationFeatureMapping',
            -rc_name           => 'default_mem',
            -analysis_capacity => 5,
            -max_retry_count  => 0,
            -flow_into => { 
                '2->A' => ['run_mapping'],
                'A->1' => ['init_parse_mapping']
            },		
        },
        {
            -logic_name => 'run_mapping',
            -module     => 'Bio::EnsEMBL::Variation::Pipeline::Remapping::RunMapping',
            -rc_name    => 'high_mem',
        },
        {
            -logic_name        => 'init_parse_mapping', 
            -module            => 'Bio::EnsEMBL::Variation::Pipeline::Remapping::InitParseMapping',
            -rc_name           => 'default_mem',
            -analysis_capacity => 5,
            -flow_into => { 
                '2->A' => ['parse_mapping'],
                'A->1' => ['init_filter_mapping']
            },    
        },
        {
            -logic_name        => 'parse_mapping',
            -module            => 'Bio::EnsEMBL::Variation::Pipeline::Remapping::ParseMapping',
            -analysis_capacity => 5,
            -rc_name           => 'default_mem',
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
            -module            => 'Bio::EnsEMBL::Variation::Pipeline::Remapping::FilterStructuralVariationFeatureMapping',
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
        {
          -logic_name => 'load_mapping',
          -module     => 'Bio::EnsEMBL::Variation::Pipeline::Remapping::LoadMapping',
          -flow_into => {
                1 => ['qc_mapping'],
            },
        },
        {
          -logic_name => 'qc_mapping',
          -module     => 'Bio::EnsEMBL::Variation::Pipeline::Remapping::StructuralVariationFeatureQC',
          -flow_into => {
                1 => ['compare_prev_assembly'],
            },
        },
        {
          -logic_name => 'compare_prev_assembly',
          -module     => 'Bio::EnsEMBL::Variation::Pipeline::Remapping::ComparePreviousAssembly',
        }
    );
   return \@analyses;
}

1;

