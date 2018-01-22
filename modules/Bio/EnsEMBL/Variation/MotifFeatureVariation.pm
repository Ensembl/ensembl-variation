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

package Bio::EnsEMBL::Variation::MotifFeatureVariation;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::MotifFeatureVariationAllele;
use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(overlap);

use base qw(Bio::EnsEMBL::Variation::RegulationVariation);

sub new {
    my $class = shift;

    my %args = @_;

    # swap a '-motif_feature' argument for a '-feature' one for the superclass
    unless($args{'-feature'} ||= delete $args{'-motif_feature'}) {
      for my $arg (keys %args) {
        if (lc($arg) eq '-motif_feature') {
          $args{'-feature'} = delete $args{$arg};
        }
      }
    }

    # call the superclass constructor
    my $self = $class->SUPER::new(%args) || return undef;
    
    # rebless the alleles from vfoas to mfvas
    map { bless $_, 'Bio::EnsEMBL::Variation::MotifFeatureVariationAllele' } 
        @{ $self->get_all_MotifFeatureVariationAlleles };
    
    return $self;
}

sub motif_feature {
    my ($self, $mf) = @_;
    return $self->SUPER::feature($mf, 'MotifFeature');
}

sub motif_feature_id {
    my $self = shift;
    unless ($self->{motif_feature_id}) {
       $self->{motif_feature_id} = $self->motif_feature->dbID;
    }
    return $self->{motif_feature_id};
}

sub regulatory_feature {
    my $self = shift;
    unless ($self->{regulatory_feature}) {
        my $funcgen_db = $self->{adaptor}->db->get_db_adaptor('funcgen');
        unless ($funcgen_db) {
            warn("Ensembl Funcgen database is not available.");
        }
        my $rfa = $funcgen_db->get_RegulatoryFeatureAdaptor;
        if ($self->motif_feature) {
            my $rf = $rfa->fetch_all_by_attribute_feature($self->motif_feature)->[0];
            $self->{regulatory_feature} = $rf;
        }
    }
    return $self->{regulatory_feature}; 
}

sub feature_stable_id {
    my $self = shift;
    unless ($self->{regulatory_feature_stable_id}) {
        $self->{regulatory_feature_stable_id} = $self->regulatory_feature->stable_id;
    }
    return $self->{regulatory_feature_stable_id};    
}

=head2 motif_name

  Arg [1]    : string $newval (optional)
               The new value to set the motif_name attribute to
  Example    : $motif_name = $obj->motif_name()
  Description: Getter/Setter for the motif_name attribute.  This is the
               name of the motif associated with this feature.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub motif_name {
  my $self = shift;
  unless ($self->{motif_name}) {
    $self->{motif_name} = $self->motif_feature->display_label;
  }
  return $self->{motif_name};
}

sub _motif_feature_seq {
    my $self = shift;
    my $mf = $self->motif_feature;
    my $mf_seq = $mf->{_variation_effect_feature_cache}->{seq} ||= $mf->seq;
    return $mf_seq;
}

sub add_MotifFeatureVariationAllele {
    my $self = shift;
    return $self->SUPER::add_VariationFeatureOverlapAllele(@_);
}

sub get_reference_MotifFeatureVariationAllele {
    my $self = shift;
    return $self->SUPER::get_reference_VariationFeatureOverlapAllele(@_);
}

sub get_all_alternate_MotifFeatureVariationAlleles {
    my $self = shift;
    return $self->SUPER::get_all_alternate_VariationFeatureOverlapAlleles(@_);
}

sub get_all_MotifFeatureVariationAlleles {
    my $self = shift;
    return $self->SUPER::get_all_VariationFeatureOverlapAlleles(@_);
}

1;
