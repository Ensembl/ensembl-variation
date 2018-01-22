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

Bio::EnsEMBL::Variation::VariationFeatureOverlap

=head1 SYNOPSIS

  use Bio::EnsEMBL::Variation::VariationFeatureOverlap;

  my $vfo = Bio::EnsEMBL::Variation::VariationFeatureOverlap->new(
    -feature      => $feature,
    -variation_feature  => $var_feat
  );

  print "consequence type: ", (join ",", @{ $vfo->consequence_type }), "\n";
  print "most severe consequence: ", $vfo->display_consequence, "\n";

=head1 DESCRIPTION

A VariationFeatureOverlap represents a VariationFeature which is in close
proximity to another Ensembl Feature. It is the superclass of feature-specific
objects such as TranscriptVariation and RegulatoryFeatureVariation, and has
methods common to all such objects. You will not normally instantiate this
class directly, instead instantiating one of the feature-specific subclasses.

=cut

package Bio::EnsEMBL::Variation::VariationFeatureOverlap;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Sequence qw(expand);
use Bio::EnsEMBL::Variation::Utils::Sequence qw(unambiguity_code);
use Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele;

use base qw(Bio::EnsEMBL::Variation::BaseVariationFeatureOverlap);

=head2 new

  Arg [-FEATURE] : 
  The Bio::EnsEMBL::Feature associated with the given VariationFeature

  Arg [-VARIATION_FEATURE] :
  The Bio::EnsEMBL::VariationFeature associated with the given Feature

  Arg [-ADAPTOR] :
  A Bio::EnsEMBL::Variation::DBSQL::VariationFeatureOverlapAdaptor

  Arg [-DISAMBIGUATE_SINGLE_NUCLEOTIDE_ALLELES] :
  A flag indiciating if ambiguous single nucleotide alleles should be disambiguated
  when constructing the VariationFeatureOverlapAllele objects, e.g. a Variationfeature
  with an allele string like 'T/M' would be treated as if it were 'T/A/C'. We limit
  ourselves to single nucleotide alleles to avoid the combinatorial explosion if we
  allowed longer alleles with potentially many ambiguous bases.

  Example : 
  my $vfo = Bio::EnsEMBL::Variation::VariationFeatureOverlap->new(
    -feature       => $feature,
    -variation_feature => $var_feat
  );

  Description: Constructs a new VariationFeatureOverlap instance given a VariationFeature
         and a Feature
  Returntype : A new Bio::EnsEMBL::Variation::VariationFeatureOverlap instance 
  Exceptions : throws unless both VARIATION_FEATURE and FEATURE are supplied, or if the
         supplied ADAPTOR is the wrong class
  Status   : Stable

=cut 

sub new {

  my $class = shift;
   
  my %args = @_;

  # swap a '-variation_feature' argument for a '-base_variation_feature' one for the superclass
  unless($args{'-base_variation_feature'} ||= delete $args{'-variation_feature'}) {
    for my $arg (keys %args) {
      if (lc($arg) eq '-variation_feature') {
        $args{'-base_variation_feature'} = delete $args{$arg};
      }
    }
  }

  my $self = $class->SUPER::new(%args);

  my (
    $adaptor, 
    $ref_feature, 
    $disambiguate_sn_alleles,
    $no_ref_check,
    $use_feature_ref,
  );

  if($Bio::EnsEMBL::Utils::Argument::NO_REARRANGE) {
    (
      $adaptor, 
      $ref_feature, 
      $disambiguate_sn_alleles,
      $no_ref_check,
      $use_feature_ref,
    ) = (
      $args{-adaptor},
      $args{-ref_feature},
      $args{-disambiguate_sn_alleles},
      $args{-no_ref_check},
      $args{-use_feature_ref},
    );
  }
  else {
    (
      $adaptor, 
      $ref_feature, 
      $disambiguate_sn_alleles,
      $no_ref_check,
      $use_feature_ref,
    ) = rearrange([qw(
      ADAPTOR 
      REF_FEATURE 
      DISAMBIGUATE_SINGLE_NUCLEOTIDE_ALLELES
      NO_REF_CHECK
      USE_FEATURE_REF
    )], %args);
  }

  my $variation_feature = $self->base_variation_feature;

  assert_ref($variation_feature, 'Bio::EnsEMBL::Variation::VariationFeature') if $Bio::EnsEMBL::Utils::Scalar::ASSERTIONS;
  assert_ref($adaptor, 'Bio::EnsEMBL::Variation::DBSQL::VariationFeatureOverlapAdaptor') if $Bio::EnsEMBL::Utils::Scalar::ASSERTIONS && $adaptor;

  $self->{adaptor} = $adaptor;

  my $ref_allele = $self->_get_ref_allele($ref_feature, $variation_feature, $no_ref_check, $use_feature_ref);

  # get raw allele hashes
  # we can do a lot of this just once per VF and cache it
  my $raw_allele_hashes = $self->_raw_allele_hashes(
    $variation_feature,
    $ref_allele,
    $disambiguate_sn_alleles,
    $use_feature_ref,
  );

  # then turn them in to VFOAs
  $self->add_VariationFeatureOverlapAllele(
    Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele->new(
      -variation_feature_overlap => $self,
      %$_
    )
  ) for @$raw_allele_hashes;
  
  return $self;
}

sub new_fast {
  my ($class, $hashref) = @_;

  # swap a variation_feature argument for a base_variation_feature one
  
  if ($hashref->{variation_feature}) {
    $hashref->{base_variation_feature} = delete $hashref->{variation_feature};
  }

  return $class->SUPER::new_fast($hashref);
}

sub dbID {
  my $self = shift;
  
  unless ($self->{dbID}) {
    # we don't really have a dbID, so concatenate all the dbIDs of our alleles

    $self->{dbID} = join '_', map { $_->dbID } @{ $self->get_all_alternate_VariationFeatureOverlapAlleles };
  }

  return $self->{dbID};
}

=head2 variation_feature

  Arg [1]  : (optional) A Bio::EnsEMBL::Variation::VariationFeature
  Description: Get/set the associated VariationFeature, lazy-loading it if required
  Returntype : Bio::EnsEMBL::Variation::VariationFeature
  Exceptions : throws if the argument is the wrong type
  Status   : Stable

=cut

sub variation_feature {
  my ($self, $variation_feature) = @_;

  if ($variation_feature) {
    assert_ref($variation_feature, 'Bio::EnsEMBL::Variation::VariationFeature') if $Bio::EnsEMBL::Utils::Scalar::ASSERTIONS;
    $self->base_variation_feature($variation_feature);
  }
  
  if (!exists($self->{base_variation_feature}) && (my $vf_id = $self->{_variation_feature_id})) {
    
    # lazy-load the VariationFeature

    if (my $adap = $self->{adaptor}) {

      if (my $vfa = $adap->db->get_VariationFeatureAdaptor) {
        ## return failed data via transcript as previous (slow otherwise especially for counts)
        my $failed = $self->adaptor->db()->include_failed_variations;
        $self->adaptor->db()->include_failed_variations(1);

        if (my $vf = $vfa->fetch_by_dbID($vf_id)) {
          $self->base_variation_feature($vf);
          delete $self->{_variation_feature_id};
        }
        ## reset fail flag
        $self->adaptor->db()->include_failed_variations($failed);
      }
    }
  }
  
  return $self->base_variation_feature;
}

sub _variation_feature_id {

  # get the dbID of the variation feature, using the VariationFeature object 
  # if we have one, or the internal hash value if we don't
  
  my $self = shift;
  
  if (my $vf = $self->{variation_feature} || $self->{base_variation_feature}) {
    return $vf->dbID;
  }
  elsif (my $id = $self->{_variation_feature_id}) {
    return $id;
  }
  else {
    return undef;
  }
}

sub get_VariationFeatureOverlapAllele_for_allele_seq {
  my ($self, $allele_seq) = @_;
  return $self->{_alleles_by_seq}->{$allele_seq};
}

=head2 add_VariationFeatureOverlapAllele

  Arg [1]  : A Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele instance
  Description: Add an allele to this VariationFeatureOverlap
  Returntype : none
  Exceptions : throws if the argument is not the expected type
  Status   : Stable

=cut

sub add_VariationFeatureOverlapAllele {
  my ($self, $vfoa) = @_;

  assert_ref($vfoa, 'Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele') if $Bio::EnsEMBL::Utils::Scalar::ASSERTIONS;

  $self->add_BaseVariationFeatureOverlapAllele($vfoa);

  $self->{_alleles_by_seq}->{ $vfoa->variation_feature_seq } = $vfoa;
}

=head2 get_reference_VariationFeatureOverlapAllele

  Description: Get the object representing the reference allele of this VariationFeatureOverlapAllele
  Returntype : Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele instance
  Exceptions : none
  Status   : Stable

=cut

sub get_reference_VariationFeatureOverlapAllele {
  my $self = shift;
  return $self->get_reference_BaseVariationFeatureOverlapAllele(@_);
}

=head2 get_all_alternate_VariationFeatureOverlapAlleles

  Description: Get a list of the alternate alleles of this VariationFeatureOverlapAllele
  Returntype : listref of Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele objects
  Exceptions : none
  Status   : Stable

=cut

sub get_all_alternate_VariationFeatureOverlapAlleles {
  my $self = shift;
  return $self->get_all_alternate_BaseVariationFeatureOverlapAlleles(@_);
}

=head2 get_all_VariationFeatureOverlapAlleles

  Description: Get a list of the all the alleles, both reference and alternate, of this
         VariationFeatureOverlap
  Returntype : listref of Bio::EnsEMBL::Variation::VariationFeatureOverlapAllele objects
  Exceptions : none
  Status   : Stable

=cut

sub get_all_VariationFeatureOverlapAlleles {
  my $self = shift;
  return $self->get_all_BaseVariationFeatureOverlapAlleles(@_);
}

sub _convert_to_sara {
  my $self = shift;
  
  my $ref_allele = $self->{reference_allele};
  $ref_allele->_convert_to_sara;
  
  $self->{alt_alleles} = [$ref_allele];
}

sub _rearrange_alleles {
  my $self = shift;
  my $keep_alleles = shift;
  
  # fix alt alleles
  my $alt_alleles = $self->{alt_alleles};
  my @new_alleles = grep {$keep_alleles->{$_->variation_feature_seq}} @$alt_alleles;
  $self->{alt_alleles} = scalar @new_alleles ? \@new_alleles : $alt_alleles;
  
  # copy to ref allele if homozygous non-ref
  $self->{reference_allele} = $self->{alt_alleles}->[0] if scalar keys %$keep_alleles == 1;
}

sub _get_ref_allele {
  my ($self, $ref_feature, $vf, $no_ref_check, $use_feature_ref) = @_;

  my $ref_allele;

  unless($no_ref_check) {
    $ref_feature ||= $vf->slice;
    $self->{ref_feature} = $ref_feature;

    my ($vf_start, $vf_end) = ($vf->start, $vf->end);

    if($vf_start > $vf_end) {
      $ref_allele = '-';
    }
    else {
      my $ss = $ref_feature->sub_Slice(
        $vf->start, 
        $vf->end, 
        $vf->strand
      );

      $ss ||= $vf->feature_Slice();

      throw("Could not fetch sub_Slice for $vf_start\-$vf_end\:".$vf->strand) unless $ss;

      $ref_allele = $ss->seq;
    }

    # throw it away if it comes back "N"
    undef $ref_allele if $ref_allele =~ /^N+$/;
  }

  return $ref_allele;
}

# this method does the raw parsing of the VFs alleles
# we cache it on VF so that we can reuse it
# as it will give the same result every time we want to create a new VFO
sub _raw_allele_hashes {
  my ($self, $vf, $ref_allele, $disambiguate_sn_alleles, $use_feature_ref) = @_;

  $vf ||= $self->base_variation_feature;

  if(!exists($vf->{_raw_allele_hashes}) || $use_feature_ref) {

    my @raw_allele_hashes;

    # get the variation feature allele string, expand it, and split it into separate alleles
    
    my $allele_string = $vf->allele_string;
    
    expand(\$allele_string);

    my @alleles = split /\//, $allele_string;

    my $vf_ref;
    if($use_feature_ref) {
      $vf_ref = $alleles[0];
      shift @alleles if defined($ref_allele);
    }
    
    $ref_allele = $alleles[0] unless defined($ref_allele);
    $ref_allele = '-' unless $ref_allele;
    
    if ($disambiguate_sn_alleles) {
      
      # if this flag is set, disambiguate any ambiguous single nucleotide alleles, so  
      # e.g. an allele string like T/M would be equivalent to an allele string of T/A/C
      # we only do this for single nucleotide alleles to avoid the combinatorial explosion
      # of long allele strings with potentially many ambiguous bases (because ensembl 
      # genomes want this functionality)

      my @possible_alleles;

      for my $allele (@alleles) {
        
        if ($allele !~ /^[ACGT-]+$/ && length($allele) == 1) {
          for my $possible ( split //, unambiguity_code($allele) ) {
            push @possible_alleles, $possible;
          }
        }
        else {
          # the allele is either unambiguous or longer than 1 nucleotide, so add it unaltered
          push @possible_alleles, $allele;
        }
      }

      @alleles = @possible_alleles;
    }

    # make sure the alleles are unique
    
    # we also want to deal with alleles like (T)0 which expand into 
    # an empty string and we want to treat this as a deletion, so 
    # we replace
    # any empty strings with '-'
    my (@new, %seen_alleles);

    foreach my $a(@alleles) {
      $a ||= '-';
      next if $seen_alleles{$a};
      push @new, $a;
      $seen_alleles{$a} = 1;
    }
    @alleles = @new;

    push @raw_allele_hashes, {
      -variation_feature_seq => $ref_allele,
      -is_reference          => 1,
      -allele_number         => 0,
      -given_ref             => $vf_ref,
    };

    # create objects representing the alternate alleles
    my $allele_number = ($vf->{_base_allele_number} || 0);

    for my $allele (@alleles) {
      
      next if $allele eq ($vf_ref || $ref_allele);

      $allele_number++;
      next if $allele eq '*' || $allele eq '<DEL:*>';

      push @raw_allele_hashes, {
        -variation_feature_seq => $allele,
        -is_reference          => 0,
        -allele_number         => $allele_number
      };
    }

    return \@raw_allele_hashes if $use_feature_ref;

    $vf->{_raw_allele_hashes} = \@raw_allele_hashes
  }

  return $vf->{_raw_allele_hashes};
}

1;
