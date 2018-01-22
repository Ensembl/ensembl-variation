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

package Bio::EnsEMBL::Variation::IntergenicVariationAllele;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::Utils::Constants qw(%OVERLAP_CONSEQUENCES);

use base qw(Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele);

sub new_fast {
    my ($self, $hashref) = @_;
    
    # swap an intergenic_variation argument for a variation_feature_overlap one

    if ($hashref->{intergenic_variation}) {
        $hashref->{variation_feature_overlap} = delete $hashref->{intergenic_variation};
    }
    
    # and call the superclass

    return $self->SUPER::new_fast($hashref);
}

sub intergenic_variation {
    my $self = shift;
    return $self->variation_feature_overlap(@_);
}

sub get_all_OverlapConsequences {
    return [ $OVERLAP_CONSEQUENCES{intergenic_variant} ];
}

1;
