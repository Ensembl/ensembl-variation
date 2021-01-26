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
package Bio::EnsEMBL::Variation::Pipeline::RegulationEffect;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::MotifFeatureVariation;
use Bio::EnsEMBL::Variation::RegulatoryFeatureVariation;
use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);

sub run {
    my $self = shift;
    my $seq_region_name = $self->param('seq_region_name');
    my $seq_region_start = $self->param('seq_region_start');
    my $seq_region_end = $self->param('seq_region_end');
    my $species = $self->param('species');
    my $disambiguate_sn_alleles = $self->param('disambiguate_single_nucleotide_alleles'); 

    my $cdba = $self->get_adaptor($species, 'core');
    $cdba->dbc->reconnect_when_lost(1);
    my $vdba = $self->get_adaptor($species, 'variation');
    $vdba->dbc->reconnect_when_lost(1);
    my $fdba = $self->get_adaptor($species, 'funcgen');
    $fdba->dbc->reconnect_when_lost(1);


    my $slice_adaptor = $cdba->get_SliceAdaptor; 
    my $regulatory_feature_adaptor = $fdba->get_RegulatoryFeatureAdaptor;
    my $motif_feature_adaptor = $fdba->get_MotifFeatureAdaptor;

    my $slice = $slice_adaptor->fetch_by_region('toplevel', $seq_region_name, $seq_region_start, $seq_region_end); 
    my @regulatory_features = grep { $_->seq_region_end <= $seq_region_end  } @{$regulatory_feature_adaptor->fetch_all_by_Slice($slice)||[]};
    $self->add_regulatory_feature_variations(\@regulatory_features);

    if ($self->param('use_experimentally_validated_mf')) {
      foreach my $rf (@regulatory_features) {
        my $motif_features = $rf->get_all_experimentally_verified_MotifFeatures();
        $self->add_motif_feature_variations($motif_features);
      }
    } else {
      my @motif_features = grep { $_->seq_region_end <= $seq_region_end  } @{$motif_feature_adaptor->fetch_all_by_Slice($slice)||[]};
      $self->add_motif_feature_variations(\@motif_features);
    }
    return;
}

sub add_regulatory_feature_variations {
  my $self = shift;
  my $regulatory_features = shift;
  my $species = $self->param('species');
  my $disambiguate_sn_alleles = $self->param('disambiguate_single_nucleotide_alleles'); 

  my $cdba = $self->get_adaptor($species, 'core');
  $cdba->dbc->reconnect_when_lost(1);
  my $vdba = $self->get_adaptor($species, 'variation');
  $vdba->dbc->reconnect_when_lost(1);


  my $rfva = $vdba->get_RegulatoryFeatureVariationAdaptor;
  $rfva->db->include_failed_variations(1);
  my $slice_adaptor = $cdba->get_SliceAdaptor; 

  foreach my $regulatory_feature (@$regulatory_features) {
    my $slice = $slice_adaptor->fetch_by_Feature($regulatory_feature) or die "Failed to get slice around RegulatoryFeature: " . $regulatory_feature->stable_id;
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
  }
}

sub add_motif_feature_variations {
  my $self = shift;
  my $motif_features = shift;

  my $species = $self->param('species');
  my $disambiguate_sn_alleles = $self->param('disambiguate_single_nucleotide_alleles'); 


  my $cdba = $self->get_adaptor($species, 'core');
  $cdba->dbc->reconnect_when_lost(1);
  my $vdba = $self->get_adaptor($species, 'variation');
  $vdba->dbc->reconnect_when_lost(1);

  my $mfva = $vdba->get_MotifFeatureVariationAdaptor;
  $mfva->db->include_failed_variations(1);
  my $slice_adaptor = $cdba->get_SliceAdaptor; 
   
  foreach my $motif_feature (@$motif_features) {
    my $slice = $slice_adaptor->fetch_by_Feature($motif_feature) or die "Failed to get slice around motif feature: " . $motif_feature->dbID;

    for my $vf ( @{ $slice->get_all_VariationFeatures }) {
      next if ($vf->allele_string !~ /^[ACGT\/]+$/i);
      my $mfv = Bio::EnsEMBL::Variation::MotifFeatureVariation->new(
        -motif_feature      => $motif_feature,
        -variation_feature  => $vf,
        -adaptor            => $mfva,
        -disambiguate_single_nucleotide_alleles => $disambiguate_sn_alleles,
      );
      if ($mfv && (scalar(@{$mfv->consequence_type}) > 0) ) {
        $mfva->store($mfv);
      }
    }
  }
}

1;
