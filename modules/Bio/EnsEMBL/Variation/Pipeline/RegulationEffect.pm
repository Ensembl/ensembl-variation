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
package Bio::EnsEMBL::Variation::Pipeline::RegulationEffect;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::MotifFeatureVariation;
use Bio::EnsEMBL::Variation::RegulatoryFeatureVariation;

use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);


sub run {
    my $self = shift;
    my $feature_id = $self->param('feature_id');
    my $feature_type = $self->param('feature_type');
    my $species = $self->param('species');

    my $disambiguate_sn_alleles = $self->param('disambiguate_single_nucleotide_alleles'); 

    my $cdba = $self->get_adaptor($species, 'core');
    my $vdba = $self->get_adaptor($species, 'variation');
    my $fdba = $self->get_adaptor($species, 'funcgen');

    my $sa = $cdba->get_SliceAdaptor; 
    my $rfa = $fdba->get_RegulatoryFeatureAdaptor;
    my $mfa = $fdba->get_MotifFeatureAdaptor;

    if ($feature_type eq 'regulatory_feature') {
        my $rfva = $vdba->get_RegulatoryFeatureVariationAdaptor;
        my $regulatory_feature = $rfa->fetch_by_stable_id($feature_id) 
            or die "Failed to fetch RegulatoryFeature for stable id: $feature_id";

        # we need to include failed variations
        $rfva->db->include_failed_variations(1);
        my $slice = $sa->fetch_by_Feature($regulatory_feature) or die "Failed to get slice around RegulatoryFeature: " . $regulatory_feature->stable_id;

        for my $vf ( @{ $slice->get_all_VariationFeatures }, @{ $slice->get_all_somatic_VariationFeatures } ) {
            my $rfv = Bio::EnsEMBL::Variation::RegulatoryFeatureVariation->new(
                -regulatory_feature => $regulatory_feature,
                -variation_feature  => $vf,
                -adaptor            => $rfva,
                -disambiguate_single_nucleotide_alleles => $disambiguate_sn_alleles,
            );

            if ($rfv && (scalar(@{$rfv->consequence_type}) > 0 )) {
                $rfva->store($rfv);
            }
        }
    } elsif ($feature_type eq 'motif_feature') {
        my $mfva = $vdba->get_MotifFeatureVariationAdaptor;
        my $motif_feature = $mfa->fetch_by_dbID($feature_id) 
            or die "Failed to fetch MotifFeature for id: $feature_id";
        my $rf = $rfa->fetch_all_by_attribute_feature($motif_feature)->[0];
        if ($rf) { 
            # we need to include failed variations
            $mfva->db->include_failed_variations(1);
            my $slice = $sa->fetch_by_Feature($motif_feature) or die "Failed to get slice around motif feature: " . $motif_feature->dbID;

            for my $vf ( @{ $slice->get_all_VariationFeatures }, @{ $slice->get_all_somatic_VariationFeatures } ) {
                my $mfv = Bio::EnsEMBL::Variation::MotifFeatureVariation->new(
                    -motif_feature      => $motif_feature,
                    -variation_feature  => $vf,
                    -adaptor            => $mfva,
                    -disambiguate_single_nucleotide_alleles => $disambiguate_sn_alleles,
                );
                if ($mfv && (scalar(@{$mfv->consequence_type}) > 0) ) {
                    $mfva->store($mfv, $rf);
                }
            }
        } else {
            $self->warning('No stable id for RF. MF: ' . $feature_id);
        }
    } elsif ($feature_type eq 'external_feature') {

    } else {    
        die "Feature type: $feature_type is not a valid argument";
    }
    return;
}



1;
