package Bio::EnsEMBL::Variation::StructuralVariationOverlapAllele;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::BaseVariationFeatureOverlapAllele);

sub new_fast {
    my ($class, $hashref) = @_;

    # swap a transcript_variation argument for a variation_feature_overlap one

    if ($hashref->{structural_variation_overlap}) {
        $hashref->{base_variation_feature_overlap} = 
            delete $hashref->{structural_variation_overlap};
    }
    
    # and call the superclass

    return $class->SUPER::new_fast($hashref);
}

1;

