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

Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele

=head1 SYNOPSIS

    use Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele;
    
    my $vfoa = Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele->new(
        -variation_feature_overlap  => $vfo,
        -variation_feature_seq      => 'A',
        -is_reference               => 0,
        -allele_number              => 1,
    );

    print "sequence with respect to the feature: ", $vfoa->feature_seq, "\n";
    print "sequence with respect to the variation feature: ", $vfoa->variation_feature_seq, "\n";
    print "consequence SO terms: ", (join ",", map { $_->SO_term } @{ $vfoa->get_all_OverlapConsequences }), "\n";

=head1 DESCRIPTION

A VariationFeatureOverlapAllele object represents a single allele of a 
VariationFeatureOverlap. It is the super-class of various feature-specific allele
classes such as TranscriptVariationAllele and RegulatoryFeatureVariationAllele and 
contains methods not specific to any particular feature type. Ordinarily you will 
not create these objects yourself, but instead you would create e.g. a 
TranscriptVariation object which will then create VariationFeatureOverlapAlleles 
based on the allele string of the associated VariationFeature. 

=cut

package Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele;

use strict;
use warnings;
use Data::Dumper;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);

use base qw(Bio::EnsEMBL::Variation::BaseVariationFeatureOverlapAllele);

our $UNAMBIGUOUS_NUCLEOTIDES = qr/^[ACGT-]+$/i;

our $ALL_NUCLEOTIDES = qr/^[ACGTUMRWSYKVHDBXN-]+$/i;

our $SPECIFIED_LENGTH = qr /(\d+) BP (INSERTION|DELETION)/i;

=head2 new

  Arg [-VARIATION_FEATURE_OVERLAP] : 
    The Bio::EnsEMBL::VariationFeatureOverlap with which this allele is 
    associated

  Arg [-VARIATION_FEATURE_SEQ] :
    The allele sequence with respect to the associated VariationFeature

  Arg [-IS_REFERENCE] :
    A flag indicating if this allele is the reference allele or not

  Arg [-ALLELE_NUMBER] :
    The order in which this allele appears in the BaseVariationFeature's
    allele string

  Example : 
    my $vfoa = Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele->new(
        -variation_feature_ovelap   => $vfo,
        -variation_feature_seq      => 'A',
        -is_reference               => 0
        -allele_number              => 1,
    );

  Description: Constructs a new VariationFeatureOverlapAllele instance given a 
               VariationFeatureOverlap and the sequence of the allele
  Returntype : A new Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele instance 
  Exceptions : throws unless both VARIATION_FEATURE_OVERLAP and VARIATION_FEATURE_SEQ 
               are supplied
  Status     : Stable

=cut 

sub new {
  my $class = shift;

  my %args = @_;

  # swap a '-variation_feature_overlap' argument for a '-base_variation_feature_overlap'
  # and a '-variation_feature' for a '-base_variation_feature' for the superclass
  unless($args{'-base_variation_feature_overlap'} ||= delete $args{'-variation_feature_overlap'}) {
    for my $arg (keys %args) {
      if (lc($arg) eq '-variation_feature_overlap') {
        $args{'-base_variation_feature_overlap'} = delete $args{$arg};
      }
    }
  }

  my $self = $class->SUPER::new(%args);

  assert_ref($self->base_variation_feature_overlap, 'Bio::EnsEMBL::Variation::VariationFeatureOverlap') if $Bio::EnsEMBL::Utils::Scalar::ASSERTIONS;

  my (
    $variation_feature_seq,
    $allele_number,
    $given_ref,
  );

  if($Bio::EnsEMBL::Utils::Argument::NO_REARRANGE) {
    (
      $variation_feature_seq,
      $allele_number,
      $given_ref,
    ) = (
      $args{-variation_feature_seq},
      $args{-allele_number},
      $args{-given_ref},
    );
  }
  else {
    (
      $variation_feature_seq,
      $allele_number,
      $given_ref,
    ) = rearrange([qw(
      VARIATION_FEATURE_SEQ
      ALLELE_NUMBER
      GIVEN_REF
    )], %args);
  }


  throw("Allele sequence required (variation "+$self->variation_feature->variation_name+")") 
    unless $variation_feature_seq;

  $self->{variation_feature_seq} = $variation_feature_seq;
  $self->{allele_number} = $allele_number;
  $self->{given_ref} = $given_ref;

  return $self;
}

sub new_fast {
    my ($class, $hashref, $strong ) = @_;
    
    # swap a variation_feature_overlap argument for a base_variation_feature_overlap one

    if ($hashref->{variation_feature_overlap}) {
        $hashref->{base_variation_feature_overlap} = delete $hashref->{variation_feature_overlap};
    }
    
    # and call the superclass

    return $class->SUPER::new_fast($hashref, $strong);
}

=head2 dbID

  Description: Get/set the dbID of this VariationFeatureOverlapAllele
  Returntype : integer
  Exceptions : none 
  Status     : Stable

=cut

sub dbID {
    my ($self, $dbID) = @_;
    $self->{dbID} = $dbID if defined $dbID;
    return $self->{dbID};
}

=head2 variation_feature_overlap

  Description: Get/set the associated VariationFeatureOverlap
  Returntype : Bio::EnsEMBL::Variation::VariationFeatureOverlap
  Exceptions : throws if the argument is the wrong type
  Status     : Stable

=cut

sub variation_feature_overlap {
    my ($self, $variation_feature_overlap) = @_;
    
    if ($variation_feature_overlap) {
        assert_ref($variation_feature_overlap, 'Bio::EnsEMBL::Variation::VariationFeatureOverlap') if $Bio::EnsEMBL::Utils::Scalar::ASSERTIONS;
    }
    
    return $self->base_variation_feature_overlap($variation_feature_overlap);
}

=head2 variation_feature

  Description: Get the associated VariationFeature
  Returntype : Bio::EnsEMBL::Variation::VariationFeature
  Exceptions : none
  Status     : Stable

=cut

sub variation_feature {
    my $self = shift;
    return $self->variation_feature_overlap->variation_feature;
}

=head2 feature_seq

  Description: Get the sequence of this allele relative to the associated Feature.
               This will be the same as the variation_feature_seq when the associated
               VariationFeature is on the same strand as the Feature, or the reverse
               complement when the strands differ.
  Returntype : string
  Exceptions : none
  Status     : At Risk

=cut

sub feature_seq {
    my $self = shift;
    
    unless ($self->{feature_seq}) {
        
        # check if we need to reverse complement the variation_feature_seq
        
        if (($self->variation_feature->strand != $self->feature->strand) && $self->seq_is_dna) {
            my $vf_seq = $self->variation_feature_seq;
            reverse_comp(\$vf_seq);
            $self->{feature_seq} = $vf_seq;
        }
        else {
            $self->{feature_seq} = $self->{variation_feature_seq};
        }
    }
    
    return $self->{feature_seq};
}

=head2 variation_feature_seq

  Args [1]   : The allele sequence relative to the VariationFeature
  Description: Get/set the sequence of this allele relative to the associated VariationFeature.
  Returntype : string
  Exceptions : none
  Status     : At Risk

=cut

sub variation_feature_seq {
    # the sequence of this allele relative to the variation feature
    my ($self, $variation_feature_seq) = @_;
    $self->{variation_feature_seq} = $variation_feature_seq if $variation_feature_seq;
    return $self->{variation_feature_seq};
}

=head2 seq_is_unambiguous_dna

  Description: identify if the sequence of this allele is unambiguous DNA 
               i.e. if we can meaningfully translate it
  Returntype : bool
  Exceptions : none
  Status     : At Risk

=cut

sub seq_is_unambiguous_dna {
    my $self = shift;

    unless (defined $self->{seq_is_unambiguous_dna}) {
        $self->{seq_is_unambiguous_dna} = 
            $self->{variation_feature_seq} =~ /$UNAMBIGUOUS_NUCLEOTIDES/ ? 1 : 0;
    }
    
    return $self->{seq_is_unambiguous_dna};
}

=head2 seq_is_dna

  Description: identify if the sequence of this allele is DNA including ambiguity 
               codes, use seq_is_unambiguous_dna to check for alleles that do not
               include ambiguity codes
  Returntype : bool
  Exceptions : none
  Status     : At Risk

=cut

sub seq_is_dna {
    my $self = shift;

    unless (defined $self->{seq_is_dna}) {
        $self->{seq_is_dna} = 
            $self->{variation_feature_seq} =~ /$ALL_NUCLEOTIDES/ ? 1 : 0;
    }
    
    return $self->{seq_is_dna};
}

=head2 seq_length

  Description: return the length of this allele sequence, this is better than
               just using length($vfoa->feature_seq) because we check if the
               sequence is valid DNA, and also look for allele strings like 
               "(3 BP INSERTION)" to determine the length
  Returntype : int or undef if we cannot determine the length
  Exceptions : none
  Status     : At Risk

=cut

sub seq_length {
    my $self = shift;

    my $seq = $self->variation_feature_seq;

    if ($self->seq_is_dna) {
        if ($seq eq '-') {
            return 0;
        }
        else {
            return length($seq);
        }
    }
    elsif ($seq =~ /$SPECIFIED_LENGTH/) {
        return $1;
    }
    
    return undef;
}

=head2 allele_string

  Description: Return a '/' delimited string of the reference allele variation_feature_seq 
               and the variation_feature_seq of this allele
  Returntype : string
  Exceptions : none
  Status     : At Risk

=cut

sub allele_string {
    my $self = shift;
    
    my $ref = $self->base_variation_feature_overlap->get_reference_VariationFeatureOverlapAllele->variation_feature_seq;
    
    # for the HGMDs and CNV probes where the alleles are artificially set to be
    # the same, just return the reference sequence
    
    if ($ref eq $self->variation_feature_seq) {
        return $ref;
    }
    else {
        return $ref.'/'.$self->variation_feature_seq;
    }
}


sub _convert_to_sara {
    my $self = shift;
    
    my $oc = Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
        'label'        => 'SARA',
        'description'  => 'Same as reference allele',
        'rank'         => '99',
        'display_term' => 'SARA',
        'SO_term'      => 'SARA',
    });
    
    $self->add_OverlapConsequence($oc);
    
    return $self;
}

# establish which parts of the given ALT sequence are different from REF
# returns listref of hashrefs
# { s => $start, e => $end }
# coords are 0-indexed
sub _get_differing_regions {
  my ($self, $bvfo) = @_;

  if(!exists($self->{_differing_regions})) {
    $bvfo ||= $self->base_variation_feature_overlap;

    my $pre = $self->_pre_consequence_predicates();
    my ($ref_length, $alt_length) = ($pre->{ref_length}, $pre->{alt_length});

    my @diff_regions;

    # obviously we don't need to do this for SNPs or anything where one of ALT or REF is 1 base or less
    # otherwise it wouldn't be a variant...
    if($alt_length > 1 && $ref_length > 1) {
      my $al1 = $bvfo->get_reference_BaseVariationFeatureOverlapAllele->variation_feature_seq;
      my $al2 = $self->variation_feature_seq;

      $al1 = '' if $al1 eq '-';
      $al2 = '' if $al2 eq '-';

      ## this code finds diffs between the sequences
      ## copied verbatim from http://www.perlmonks.org/?node_id=882593
      my ($m, @pos);
      $m = $al1 ^ $al2;
      push @pos, pos( $m ) while $m =~ m{(?=([^\x00]))}g;

      # now we join them into contiguous regions
      my ($region, $prev);

      foreach my $pos(@pos) {
        if($region && $pos == $prev + 1) {
          $region->{e} = $pos;
        }
        else {
          push @diff_regions, $region if $region;

          $region = {
            s => $pos,
            e => $pos
          };
        }

        $prev = $pos;
      }

      push @diff_regions, $region if $region;
    }
    else {
      push @diff_regions, { s => 0, e => $ref_length - 1 };
    }

    $self->{_differing_regions} = \@diff_regions;
  }

  return $self->{_differing_regions};
}

1;