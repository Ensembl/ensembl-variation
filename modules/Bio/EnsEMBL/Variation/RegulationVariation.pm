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

package Bio::EnsEMBL::Variation::RegulationVariation;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::VariationFeatureOverlap);

sub feature_label {
    my ($self, $feature_label) = @_;
    $self->{feature_label} = $feature_label if $feature_label;
    return $self->{feature_label};
}

sub target_feature {
    # XXX: fetch the target feature
}

sub target_feature_stable_id {
    my ($self, $target_feature_stable_id) = @_;
    $self->{target_feature_stable_id} = $target_feature_stable_id if $target_feature_stable_id;
    return $self->{target_feature_stable_id};
}

1;
