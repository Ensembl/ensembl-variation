=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2024] EMBL-European Bioinformatics Institute

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

package Bio::EnsEMBL::Variation::StructuralVariationOverlapAllele;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);

use base qw(Bio::EnsEMBL::Variation::BaseVariationFeatureOverlapAllele);

=head2 new

  Arg [-STRUCTURAL_VARIATION_OVERLAP] : 
    The Bio::EnsEMBL::StructuralVariationOverlap with which this allele is 
    associated

  Arg [-SYMBOLIC_ALLELE] :
    The symbolic allele string

  Arg [-BREAKEND] :
    Breakend information 

  Example : 
    my $vfoa = Bio::EnsEMBL::Variation::StructuralVariationOverlapAllele->new(
        -structural_variation_overlap   => $svfo,
        -symbolic_allele                => 'N[8:56445865[',
        -breakend                       => $breakend
    );

  Description: Constructs a new StructuralVariationOverlapAllele instance given a 
               StructuralVariationOverlap and the sequence of the allele
  Returntype : A new Bio::EnsEMBL::Variation::StructuralVariationOverlapAllele instance 
  Exceptions : throws unless STRUCTURAL_VARIATION_OVERLAP is supplied
  Status     : Stable

=cut 

sub new {
    my $class = shift;

    my %args = @_;

    # swap a '-structural_variation_overlap' argument for a '-base_variation_feature_overlap'
    # and a '-variation_feature' for a '-base_variation_feature' for the superclass
    unless($args{'-base_variation_feature_overlap'} ||= delete $args{'-structural_variation_overlap'}) {
        for my $arg (keys %args) {
            if (lc($arg) eq '-structural_variation_overlap') {
                $args{'-base_variation_feature_overlap'} = delete $args{$arg};
            }
        }
    }

  my $self = $class->SUPER::new(%args);

  assert_ref($self->base_variation_feature_overlap, 'Bio::EnsEMBL::Variation::StructuralVariationOverlap') if $Bio::EnsEMBL::Utils::Scalar::ASSERTIONS;

  my (
    $symbolic_allele,
    $allele_number,
    $is_reference,
    $breakend,
  );

  if($Bio::EnsEMBL::Utils::Argument::NO_REARRANGE) {
    (
      $symbolic_allele,
      $breakend,
    ) = (
      $args{-symbolic_allele},
      $args{-breakend},
    );
  }
  else {
    (
      $symbolic_allele,
      $breakend,
    ) = rearrange([qw(
      SYMBOLIC_ALLELE
      BREAKEND
    )], %args);
  }

  $self->{symbolic_allele} = $symbolic_allele;
  $self->{breakend} = $breakend if defined $breakend;

  return $self;
}

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

=head2 structural_variation_overlap

  Description: Get the associated StructuralVariationOverlap
  Returntype : Bio::EnsEMBL::Variation::StructuralVariationOverlap
  Exceptions : none
  Status     : Stable

=cut

sub structural_variation_overlap {
    my ($self, $svo) = @_;
    if ($svo) {
        assert_ref($svo, 'Bio::EnsEMBL::Variation::StructuralVariationOverlap');
    }
    return $self->base_variation_feature_overlap($svo);
}


=head2 structural_variation_feature

  Description: Get the associated StructuralVariationFeature
  Returntype : Bio::EnsEMBL::Variation::StructuralVariationFeature
  Exceptions : none
  Status     : Stable

=cut

sub structural_variation_feature {
    my $self = shift;
    return $self->structural_variation_overlap->structural_variation_feature;
}

=head2 symbolic_allele

  Args [1]   : The symbolic allele string
  Description: Get/set the string of this symbolic allele.
  Returntype : string
  Exceptions : none
  Status     : At Risk

=cut

sub symbolic_allele {
    my ($self, $sva) = @_;
    $self->{symbolic_allele} = $sva if $sva;
    return $self->{symbolic_allele};
}

=head2 breakend

  Args [1]   : The breakend information for this allele
  Description: Get/set the breakend for this allele
  Returntype : string
  Exceptions : none
  Status     : At Risk

=cut

sub breakend {
    my ($self, $breakend) = @_;
    $self->{breakend} = $breakend if $breakend;
    return $self->{breakend};
}

# required for intron overlap checking
sub _get_differing_regions {
  my $self = shift;
  my $bvf = $self->base_variation_feature;
  return $self->{_differing_regions} ||= [{ s => 0, e => ($bvf->{end} - $bvf->{start}) }];
}

1;

