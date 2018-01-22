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

=head1 NAME

Bio::EnsEMBL::Variation::BaseTranscriptVariationAllele

=head1 SYNOPSIS

    use Bio::EnsEMBL::Variation::BaseTranscriptVariationAllele;

=head1 DESCRIPTION

An helper class for representing an overlap of a Transcript and a
Variation allele (either sequence or structural). Should not be invoked
directly.

=cut

package Bio::EnsEMBL::Variation::BaseTranscriptVariationAllele;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(overlap _intron_overlap);

use base qw(Bio::EnsEMBL::Variation::BaseVariationFeatureOverlapAllele);

=head2 base_transcript_variation

  Description: Get/set the associated BaseTranscriptVariation
  Returntype : Bio::EnsEMBL::Variation::BaseTranscriptVariation
  Exceptions : throws if the argument is the wrong type
  Status     : Stable

=cut

sub base_transcript_variation {
    my ($self, $btv) = @_;
    assert_ref($btv, 'Bio::EnsEMBL::Variation::BaseTranscriptVariation') if $btv;
    return $self->base_variation_feature_overlap($btv);
}

=head2 transcript

  Description: Get the associated Transcript
  Returntype : Bio::EnsEMBL::Transcript
  Exceptions : none
  Status     : Stable

=cut

sub transcript {
    my $self = shift;
    return $self->base_transcript_variation->transcript;
}

=head2 base_variation_feature

  Description: Get the associated BaseVariationFeature
  Returntype : Bio::EnsEMBL::Variation::BaseVariationFeature
  Exceptions : none
  Status     : Stable

=cut

sub base_variation_feature {
    my $self = shift;
    return $self->base_transcript_variation->base_variation_feature;
}

sub _intron_effects {
  my ($self, $feat, $tv, $vf) = @_;

  # internal method used by Bio::EnsEMBL::Variation::Utils::VariationEffect
  # when calculating various consequence types
  
  # this method is a major bottle neck in the effect calculation code so 
  # we cache results and use local variables instead of method calls where
  # possible to speed things up - caveat bug-fixer!
  
  unless ($self->{_intron_effects}) {
    
    my $intron_effects = {};

    $tv ||= $self->base_variation_feature_overlap;
    $vf ||= $self->base_variation_feature;
    my $vf_start = $vf->{start};
    
    foreach my $region(@{$self->_get_differing_regions($tv)}) {
      my ($r_start, $r_end) = ($vf_start + $region->{s}, $vf_start + $region->{e});

      my $insertion = $r_start == $r_end + 1;
      
      my ($min, $max) = $r_start > $r_end ? ($r_end, $r_start) : ($r_start, $r_end);

      # check introns themselves
      for my $intron (@{$tv->_overlapped_introns($min, $max)}) {

        my $intron_start = $intron->{start};
        my $intron_end   = $intron->{end};

        # check frameshift intron
        # this hash key gets set when the trees are created
        if($intron->{_frameshift} && overlap($r_start, $r_end, $intron_start, $intron_end)) {
          $intron_effects->{within_frameshift_intron} = 1;
          next;
        }

        if (
          overlap($r_start, $r_end, $intron_start+2, $intron_end-2) or 
          ($insertion && ($r_start == $intron_start+2 || $r_end == $intron_end-2))
        ) {
          $intron_effects->{intronic} = 1;
        }
      }

      # now check intron boundaries
      for my $intron (@{$tv->_overlapped_introns_boundary($min, $max)}) {

        my $intron_start = $intron->{start};
        my $intron_end   = $intron->{end};

        # check frameshift intron
        if($intron->{_frameshift} && overlap($r_start, $r_end, $intron_start, $intron_end)) {
          $intron_effects->{within_frameshift_intron} = 1;
          next;
        }

        if (overlap($r_start, $r_end, $intron_start, $intron_start+1)) {
          $intron_effects->{start_splice_site} = 1;
        }
        
        if (overlap($r_start, $r_end, $intron_end-1, $intron_end)) {
          $intron_effects->{end_splice_site} = 1;
        }
        
        # the definition of splice_region (SO:0001630) is "within 1-3 bases 
        # of the exon or 3-8 bases of the intron", the intron start is the 
        # first base of the intron so we only need to add or subtract 7 from 
        # it to get the correct coordinate. We also need to special case 
        # insertions between the edge of an exon and a donor or acceptor site
        # and between a donor or acceptor site and the intron            
        $intron_effects->{splice_region} = _intron_overlap($r_start, $r_end, $intron_start, $intron_end, $insertion)
          unless $intron_effects->{start_splice_site} or $intron_effects->{end_splice_site};
      }
    }
      
    $self->{_intron_effects} = $intron_effects;       
  }

  return $self->{_intron_effects};
}

1;