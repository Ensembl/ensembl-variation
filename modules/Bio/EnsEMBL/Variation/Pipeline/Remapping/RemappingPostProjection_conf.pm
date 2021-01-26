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

package Bio::EnsEMBL::Variation::Pipeline::Remapping::RemappingPostProjection_conf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Variation::Pipeline::Remapping::RemappingVariationFeature_conf');

sub default_options {
  my ($self) = @_;

# The hash returned from this function is used to configure the
# pipeline, you can supply any of these options on the command
# line to override these default values.

# You shouldn't need to edit anything in this file other than
# these values, if you find you do need to then we should probably
# make it an option here, contact the variation team to discuss
# this - patches are welcome!
#
  return {
    %{ $self->SUPER::default_options() },   # inherit from parent
    mode                    => 'remap_variation_feature',
    debug                   => 0,
    debug_sequence_name     => 3,
    flank_seq_length        => 150,
    algn_score_threshold    => 0.95,
    max_map_weight          => 5,
    use_prior_for_filtering => 1,
    map_to_chrom_only       => 1,
    entries_per_file        => 200000,
  };
}

sub pipeline_wide_parameters {
    my ($self) = @_;
    return {
        %{$self->SUPER::pipeline_wide_parameters},
        flank_seq_length             => $self->o('flank_seq_length'),	
        feature_table                => $self->o('feature_table'),
        algn_score_threshold         => $self->o('algn_score_threshold'),
        max_map_weight               => $self->o('max_map_weight'),
        use_prior_for_filtering      => $self->o('use_prior_for_filtering'),
        map_to_chrom_only            => $self->o('map_to_chrom_only'),
        entries_per_file             => $self->o('entries_per_file'),
        debug                        => $self->o('debug'),
        debug_sequence_name          => $self->o('debug_sequence_name'),
    };
}

sub resource_classes {
  my ($self) = @_;
  return $self->SUPER::resource_classes;
}

sub pipeline_analyses {
  my ($self) = @_;
  return $self->SUPER::pipeline_analyses;
}

1;

