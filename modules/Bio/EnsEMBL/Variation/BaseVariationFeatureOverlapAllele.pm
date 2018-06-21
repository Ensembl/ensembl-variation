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

Bio::EnsEMBL::Variation::BaseVariationFeatureOverlapAllele

=head1 SYNOPSIS

    use Bio::EnsEMBL::Variation::BaseVariationFeatureOverlapAllele;
    
    my $bvfoa = Bio::EnsEMBL::Variation::BaseVariationFeatureOverlapAllele->new(
        -base_variation_feature_overlap => $bvfo,
        -is_reference                   => 0,
    );

    print "consequence SO terms: ", (join ",", map { $_->SO_term } @{ $bvfoa->get_all_OverlapConsequences }), "\n";

=head1 DESCRIPTION

A BaseVariationFeatureOverlapAllele object represents a single allele of a 
BaseVariationFeatureOverlap. It is the super-class of variation feature specific
classes such as VariationFeatureOverlapAllele and StructuralVariationOverlapAllele 
and contains methods not specific to any particular variation feature type. 
Ordinarily you will not create these objects yourself, but instead you would 
create one of the more specific subclasses.

=cut

package Bio::EnsEMBL::Variation::BaseVariationFeatureOverlapAllele;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::EnsEMBL::Variation::Utils::Constants qw(%OVERLAP_CONSEQUENCES $DEFAULT_OVERLAP_CONSEQUENCE);
use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(overlap);
use Scalar::Util qw(weaken);

our @SORTED_OVERLAP_CONSEQUENCES = sort {$a->tier <=> $b->tier} values %OVERLAP_CONSEQUENCES;

# get which keys we care about from OCs
our %OC_INCLUDE_KEYS = map {$_ => 1} map {keys %{$_->{include} || {}}} @SORTED_OVERLAP_CONSEQUENCES;
$OC_INCLUDE_KEYS{feature_class} = 1;
$OC_INCLUDE_KEYS{vf_class} = 1;

=head2 new

  Arg [-BASE_VARIATION_FEATURE_OVERLAP] : 
    The Bio::EnsEMBL::BaseVariationFeatureOverlap with which this allele is 
    associated

  Arg [-IS_REFERENCE] :
    A flag indicating if this allele is the reference allele or not

  Example : 
    my $bvfoa = Bio::EnsEMBL::Variation::BaseVariationFeatureOverlapAllele->new(
        -base_variation_feature_overlap  => $bvfo,
        -is_reference                   => 0
    );

  Description: Constructs a new BaseVariationFeatureOverlapAllele instance given a 
               BaseVariationFeatureOverlap and a flag indicating if this is the 
               reference allele
  Returntype : A new Bio::EnsEMBL::Variation::BaseVariationFeatureOverlapAllele instance 
  Exceptions : throws unlessBASE_VARIATION_FEATURE_OVERLAP is supplied
  Status     : Stable

=cut 

sub new {
  my $class = shift;

  my (
    $base_variation_feature_overlap,
    $is_reference
  );

  if($Bio::EnsEMBL::Utils::Argument::NO_REARRANGE) {
    my %args = @_;
    (
      $base_variation_feature_overlap,
      $is_reference
    ) = (
      $args{-base_variation_feature_overlap},
      $args{-is_reference},
    );
  }
  else {
    (
      $base_variation_feature_overlap,
      $is_reference
    ) = rearrange([qw(
      BASE_VARIATION_FEATURE_OVERLAP
      IS_REFERENCE
    )], @_);
  }

  assert_ref($base_variation_feature_overlap, 'Bio::EnsEMBL::Variation::BaseVariationFeatureOverlap') if $Bio::EnsEMBL::Utils::Scalar::ASSERTIONS;
  
  my $self = bless {
    base_variation_feature_overlap  => $base_variation_feature_overlap,
    is_reference                    => $is_reference,
    }, $class;

  # avoid a memory leak, because the bvfo also has a reference to us
  weaken $self->{base_variation_feature_overlap};

  return $self;
}

sub new_fast {
    my ($class, $hashref, $strong) = @_;
    my $self = bless $hashref, $class;
    # avoid a memory leak, because the bvfo also has a reference to us
    weaken $self->{base_variation_feature_overlap} if $self->{base_variation_feature_overlap} && !defined $strong;
    return $self;
}

=head2 base_variation_feature_overlap

  Description: Get/set the associated BaseVariationFeatureOverlap
  Returntype : Bio::EnsEMBL::Variation::BaseVariationFeatureOverlap
  Exceptions : throws if the argument is the wrong type
  Status     : Stable

=cut

sub base_variation_feature_overlap {
    my ($self, $bvfo) = @_;

    if ($bvfo) {
        assert_ref($bvfo, 'Bio::EnsEMBL::Variation::BaseVariationFeatureOverlap') if $Bio::EnsEMBL::Utils::Scalar::ASSERTIONS;
        $self->{base_variation_feature_overlap} = $bvfo;
        # avoid a memory leak, because the bvfo also has a reference to us
        weaken $self->{base_variation_feature_overlap};
    }

    return $self->{base_variation_feature_overlap};
}

=head2 base_variation_feature

  Description: Get the associated BaseVariationFeature
  Returntype : Bio::EnsEMBL::Variation::BaseVariationFeature
  Exceptions : none
  Status     : Stable

=cut

sub base_variation_feature {
    my $self = shift;
    return $self->base_variation_feature_overlap->base_variation_feature(@_);
}

=head2 feature

  Description: Get the associated Feature
  Returntype : Bio::EnsEMBL::Feature (or relevant subclass)
  Exceptions : none
  Status     : Stable

=cut

sub feature {
    my $self = shift;
    return $self->base_variation_feature_overlap->feature(@_);
}

=head2 is_reference

  Args [1]   : A boolean value 
  Description: Get/set a flag indicating if this allele is the reference allele
  Returntype : bool
  Exceptions : none
  Status     : Stable

=cut

sub is_reference {
    my ($self, $is_reference) = @_;
    $self->{is_reference} = $is_reference if defined $is_reference;
    return $self->{is_reference};
}

=head2 allele_number

  Description: Getter/setter for the allele number; represents the order of this allele
               in the BaseVariationFeature's allele string
  Returntype : int
  Exceptions : none
  Status     : At Risk

=cut

sub allele_number {
  my $self = shift;

  $self->{allele_number} = shift if @_;

  return $self->{allele_number};
}

=head2 get_all_OverlapConsequences

  Description: Get a list of all the OverlapConsequences of this allele, calculating them 
               on the fly if necessary
  Returntype : listref of Bio::EnsEMBL::Variation::OverlapConsequence objects
  Exceptions : none
  Status     : Stable

=cut

sub get_all_OverlapConsequences {
  my $self = shift;
  
  unless ($self->{overlap_consequences}) {

    # calculate consequences on the fly
    
    my $cons = [];
    
    my $assigned_tier;
    
    # get references to frequently used objects
    # by passing these to the predicates here we reduce the number of method calls dramatically
    my $bvfo = $self->base_variation_feature_overlap;
    my $bvf  = $bvfo->base_variation_feature;
    my $feat = $bvfo->feature;
    
    my $pre = $self->_pre_consequence_predicates($feat, $bvfo, $bvf);

    # loop over all the consequences determined for this pre hash
    OC: for my $oc (@{$self->_get_oc_list($pre)}) {

      last if $assigned_tier && $oc->{tier} > $assigned_tier;
      
      if($oc->predicate->($self, $feat, $bvfo, $bvf)) {
        push @$cons, $oc;

        my $tier = $oc->{tier};
        if($tier <= 2) {
          $assigned_tier = $tier;
        }
      }
    }

    $cons = [$DEFAULT_OVERLAP_CONSEQUENCE] unless @$cons;

    $self->{overlap_consequences} = $cons;
  }
  
  return $self->{overlap_consequences};
}

sub _get_oc_list {
  my ($self, $pre) = @_;

  my $cache = $main::_VEP_CACHE->{oc} ||= {};
  my $digest = $pre->{_digest};

  unless(exists($cache->{$digest})) {
    my $dummy_feat = bless {}, $pre->{feature_class};
    my $dummy_vf   = bless {}, $pre->{vf_class};

    $cache->{$digest} ||= [
      map {$_->{SO_term}}
      grep {
        $_->{predicate} &&
        !$self->_skip_oc($_, $pre) &&
        $dummy_vf->isa($_->{variant_feature_class}) &&
        $dummy_feat->isa($_->{feature_class})
      }
      @SORTED_OVERLAP_CONSEQUENCES
    ];
  }

  return [map {$OVERLAP_CONSEQUENCES{$_}} @{$cache->{$digest}}];
}

sub _skip_oc {
  my ($self, $oc, $pre) = @_;

  # check include to hash to determine whether to run this predicate
  if(my $inc = $oc->{include}) {
    foreach my $key(keys %$inc) {
      if($inc->{$key} == 0) {
        return 1 if exists($pre->{$key}) && $pre->{$key} ne $inc->{$key};
      }
      else {
        return 1 unless exists($pre->{$key}) && $pre->{$key} eq $inc->{$key};
      }
    }
  }

  return 0;
}

# pre-calculate a bunch of attributes
# used to work out whether to execute predicates in get_all_OverlapConsequences
sub _pre_consequence_predicates {
  my ($self, $feat, $bvfo, $bvf) = @_;
  
  unless(exists($self->{pre_consequence_predicates})) {
    $bvfo ||= $self->base_variation_feature_overlap;
    $bvf  ||= $bvfo->base_variation_feature;
    $feat ||= $bvfo->feature;
    
    my $preds = {};

    # this string will uniquely identify the combination of keys and values in preds
    # we can use it to cache a list of predicate methods we should run on similar variants
    # so it doesn't have to be calculated every time
    my $pred_digest = '';
    
    # this _update_preds method stores the key=value pair on $preds
    # and also appends these as string to $pred_digest
    $self->_update_preds($preds, 'feature_class', ref($feat), \$pred_digest);
    $self->_update_preds($preds, 'vf_class', ref($bvf), \$pred_digest);
    
    # get BVF (variant genomic location) preds, copy to "main"
    $bvf->{pre_consequence_predicates} ||= $self->_bvf_preds($feat, $bvfo, $bvf, $preds);
    @$preds{keys %{$bvf->{pre_consequence_predicates}}} = values %{$bvf->{pre_consequence_predicates}};
    $pred_digest .= $bvf->{pre_consequence_predicates}->{_digest};
    
    # get BVFO (variant-feature location) preds, copy to "main"
    $bvfo->{pre_consequence_predicates} ||= $self->_bvfo_preds($feat, $bvfo, $bvf, $preds);
    @$preds{keys %{$bvfo->{pre_consequence_predicates}}} = values %{$bvfo->{pre_consequence_predicates}};
    $pred_digest .= $bvfo->{pre_consequence_predicates}->{_digest};
    
    ## allele-specific type for non-SVs
    unless($preds->{sv}) {
      my $vf_seq = $self->variation_feature_seq();
      $vf_seq = '' if $vf_seq eq '-';
      my $alt_length = length($vf_seq);
      $self->_update_preds($preds, 'alt_length', $alt_length, \$pred_digest);
      
      my $ref_length = $preds->{ref_length};
    
      # most variants we see will be SNPs, so test this first
      if($alt_length == $ref_length) {
        $self->_update_preds($preds, 'snp', 1, \$pred_digest);
      }
      elsif( $ref_length > $alt_length ) {
        $self->_update_preds($preds, 'deletion', 1, \$pred_digest);
        $self->_update_preds($preds, 'decrease_length', 1, \$pred_digest);
      }
      elsif( $ref_length < $alt_length ) {
        $self->_update_preds($preds, 'insertion', 1, \$pred_digest);
        $self->_update_preds($preds, 'increase_length', 1, \$pred_digest);
      }
    }
      
    # log the digest
    $preds->{_digest} = $pred_digest;

    # copy to self
    $self->{pre_consequence_predicates} = $preds;
  }
  
  return $self->{pre_consequence_predicates};
}

sub _bvf_preds {
  my ($self, $feat, $bvfo, $bvf, $preds) = @_;

  my $bvf_preds = {};
  my $pred_digest = '';

  my ($vf_start, $vf_end) = ($bvf->{start}, $bvf->{end});

  $self->_update_preds($bvf_preds, 'complete_overlap', 1, \$pred_digest)
    if $vf_start <= $feat->{start} && $vf_end >= $feat->{end};
  
  my $is_sv = $bvf->isa('Bio::EnsEMBL::Variation::StructuralVariationFeature') ? 1 : 0;
  $self->_update_preds($bvf_preds, 'sv', $is_sv, \$pred_digest);
  
  # use SO term to determine class
  if($is_sv) {
    my $class_SO_term = $bvf->class_SO_term;

    if($class_SO_term =~ /deletion/) {
      $self->_update_preds($bvf_preds, 'deletion', 1, \$pred_digest);
      $self->_update_preds($bvf_preds, 'decrease_length', 1, \$pred_digest);
    }
    elsif($class_SO_term =~ /insertion|duplication/) {
      $self->_update_preds($bvf_preds, 'insertion', 1, \$pred_digest);
      $self->_update_preds($bvf_preds, 'increase_length', 1, \$pred_digest);
    }
  }
  
  # otherwise for sequence variants, log the reference length here
  # we'll use it later to determine class per-allele
  else {
    $self->_update_preds($bvf_preds, 'ref_length', ($vf_end - $vf_start) + 1, \$pred_digest);
  }

  $bvf_preds->{_digest} = $pred_digest;
  
  return $bvf_preds;
}

sub _bvfo_preds {
  my ($self, $feat, $bvfo, $bvf, $preds) = @_;

  my $bvfo_preds = {};
  my $pred_digest = '';

  my ($vf_start, $vf_end) = ($bvf->{start}, $bvf->{end});  
  my ($min_vf, $max_vf) = $vf_start > $vf_end ? ($vf_end, $vf_start) : ($vf_start, $vf_end);
  
  # within feature
  my $wf = overlap($vf_start, $vf_end, $feat->{start}, $feat->{end}) ? 1 : 0;
  $self->_update_preds($bvfo_preds, 'within_feature', $wf, \$pred_digest);
  
  # use a complex if/else structure to avoid executing unnecessary code
  if($bvfo_preds->{within_feature}) {
    
    my ($tr_coding_start, $tr_coding_end);
    
    if($preds->{feature_class} eq 'Bio::EnsEMBL::Transcript') {
      
      # get coding region coords
      ($tr_coding_start, $tr_coding_end) = ($feat->coding_region_start, $feat->coding_region_end);
      
      # store biotype
      $self->_update_preds($bvfo_preds, $feat->biotype, 1, \$pred_digest);

      # find any overlapping introns, intron boundary regions and exons
      $self->_update_preds(
        $bvfo_preds,
        'intron',
        scalar @{$bvfo->_overlapped_introns($min_vf, $max_vf)} ? 1 : 0,
        \$pred_digest
      );

      $self->_update_preds(
        $bvfo_preds,
        'intron_boundary',
        scalar @{$bvfo->_overlapped_introns_boundary($min_vf, $max_vf)} ? 1 : 0,
        \$pred_digest
      );

      $self->_update_preds(
        $bvfo_preds,
        'exon',
        scalar @{$bvfo->_overlapped_exons($min_vf, $max_vf)} ? 1 : 0,
        \$pred_digest
      );
    }
    
    # coding transcript specific
    if($tr_coding_start && $tr_coding_end) {

      # completely within UTR
      if(!overlap($min_vf, $max_vf, $tr_coding_start, $tr_coding_end)) {
        $self->_update_preds($bvfo_preds, 'utr', 1, \$pred_digest);
      }
    
      else {
    
        # partially within UTR
        if($min_vf < $tr_coding_start || $max_vf > $tr_coding_end) {
          $self->_update_preds($bvfo_preds, 'utr', 1, \$pred_digest);
        }

        # shortcut if trees used
        if(exists($bvfo_preds->{exon}) && !$bvfo_preds->{exon}) {
          $self->_update_preds($bvfo_preds, 'non_coding', 1, \$pred_digest);
        }

        else {
          my $cds_coords = $bvfo->cds_coords;
          
          # compeletely within coding region
          if(scalar @$cds_coords == 1 && $cds_coords->[0]->isa('Bio::EnsEMBL::Mapper::Coordinate')) {
            $self->_update_preds($bvfo_preds, 'coding', 1, \$pred_digest);
          }
          
          else {
            
            # completely in non-coding region
            if(scalar @$cds_coords == 0) {
              $self->_update_preds($bvfo_preds, 'non_coding', 1, \$pred_digest);
            }
            
            # partly in coding region
            else {
              $self->_update_preds($bvfo_preds, 'coding', 1, \$pred_digest);
            }
          } 
        }
      }
    }
    
    elsif($preds->{feature_class} eq 'Bio::EnsEMBL::Transcript') {
      $self->_update_preds($bvfo_preds, 'non_coding', 1, \$pred_digest);
    }
  }

  $bvfo_preds->{_digest} = $pred_digest;

  return $bvfo_preds;
}

sub _update_preds {
  my ($self, $preds, $key, $value, $pred_digest) = @_;
  $preds->{$key} = $value;
  $$pred_digest .= $key.$value if $OC_INCLUDE_KEYS{$key};
}

=head2 add_OverlapConsequence

  Arg [1]    : Bio::EnsEMBL::Variation::OverlapConsequence instance
  Description: Add an OverlapConsequence to this allele's list 
  Returntype : none
  Exceptions : throws if the argument is the wrong type
  Status     : Stable

=cut

sub add_OverlapConsequence {
    my ($self, $oc) = @_;
    assert_ref($oc, 'Bio::EnsEMBL::Variation::OverlapConsequence') if $Bio::EnsEMBL::Utils::Scalar::ASSERTIONS;
    $self->{overlap_consequences} ||= [];
    push @{ $self->{overlap_consequences} }, $oc;
}

sub SO_isa {
    my ($self, $query) = @_;
    
    if (my $adap = $self->base_variation_feature_overlap->{adaptor}) {
        if (my $ota = $adap->db->dnadb->get_OntologyTermAdaptor) {
            my $term = $ota->fetch_by_accession();
            my @parents = $ota->fetch_by_child_term($term);
        }
    }
    
    for my $cons (@{ $self->get_all_OverlapConsequences }) {
        if ($cons->SO_term eq $query) {
            return 1;
        }
    } 
}

1;

