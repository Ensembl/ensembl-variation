=head1 LICENSE

 Copyright (c) 1999-2013 The European Bioinformatics Institute and
 Genome Research Limited.  All rights reserved.

 This software is distributed under a modified Apache license.
 For license details, please see

   http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <dev@ensembl.org>.

 Questions may also be sent to the Ensembl help desk at
 <helpdesk@ensembl.org>.

=cut

package Bio::EnsEMBL::Variation::IntergenicStructuralVariationAllele;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::Utils::Constants qw(%OVERLAP_CONSEQUENCES);

use base qw(Bio::EnsEMBL::Variation::StructuralVariationOverlapAllele);

sub new_fast {
    my ($self, $hashref) = @_;
    
    # swap an intergenic_structural_variation argument for a variation_feature_overlap one

    if ($hashref->{intergenic_structural_variation}) {
        $hashref->{variation_feature_overlap} = delete $hashref->{intergenic_structural_variation};
    }
    
    # and call the superclass

    return $self->SUPER::new_fast($hashref);
}

sub intergenic_structural_variation {
    my $self = shift;
    return $self->variation_feature_overlap(@_);
}

sub get_all_OverlapConsequences {
    return [ $OVERLAP_CONSEQUENCES{intergenic_variant} ];
}

1;
