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

package Bio::EnsEMBL::Variation::Pipeline::Remapping::Remapping_conf;

use strict;
use warnings;

use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
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
        },    
        hive_force_init         => 1,
        hive_use_param_stack    => 0,
        hive_use_triggers       => 0,
        hive_auto_rebalance_semaphores => 0,  # do not attempt to rebalance semaphores periodically by default
        hive_no_init            => 0, # setting it to 1 will skip pipeline_create_commands (useful for topping up)
        hive_root_dir           => $ENV{'HOME'} . '/bin/ensembl-hive',
        ensembl_cvs_root_dir    => $ENV{'HOME'} . '/bin',
        debug                   => 0,
        run_variant_qc          => 1,
        use_fasta_files         => 0,
        flank_seq_length        => 150,
        algn_score_threshold    => 0.95,
        max_map_weight          => 5,
        use_prior_for_filtering => 1,
        map_to_chrom_only       => 1,
        entries_per_file        => 200000,
        dump_multi_map          => 1,
        skip_patch_comparison   => 0,

        bwa                     => 'bwa',
        samtools                => 'samtools',

        registry_file_oldasm    => $self->o('pipeline_dir') . '/ensembl.registry.oldasm',
        registry_file_newasm    => $self->o('pipeline_dir') . '/ensembl.registry.newasm',

        pipeline_dir            => $self->o('pipeline_dir'),
        old_assembly            => $self->o('pipeline_dir') . '/old_assembly',
        new_assembly            => $self->o('pipeline_dir') . '/new_assembly',

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
            -dbname => $self->o('hive_db_name'),
            -driver => 'mysql',
        },
    };
}

sub pipeline_wide_parameters {
  my ($self) = @_;
    return {
      %{$self->SUPER::pipeline_wide_parameters},          # here we inherit anything from the base class
        registry_file_oldasm         => $self->o('registry_file_oldasm'),
        registry_file_newasm         => $self->o('registry_file_newasm'),
        dump_features_dir            => $self->o('dump_features'),
        dump_mapped_features_dir     => $self->o('dump_mapped_features'),
        filtered_mappings_dir        => $self->o('filtered_mappings'),
        load_features_dir            => $self->o('load_features'),
        statistics_dir               => $self->o('statistics'),
        fasta_files_dir              => $self->o('fasta_files'),
        bam_files_dir                => $self->o('bam_files'),
        mapping_results_dir          => $self->o('mapping_results'),
        old_assembly_dir             => $self->o('old_assembly'),
        new_assembly_dir             => $self->o('new_assembly'),
        species                      => $self->o('species'),
        bwa                          => $self->o('bwa'),
        samtools                     => $self->o('samtools'),
        pipeline_dir                 => $self->o('pipeline_dir'),
        flank_seq_length             => $self->o('flank_seq_length'),	
        use_fasta_files              => $self->o('use_fasta_files'),
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
            'extra_mem'   => { 'LSF' => '-R"select[mem>12500] rusage[mem=12500]" -M12500'},			
    };
}


1;

