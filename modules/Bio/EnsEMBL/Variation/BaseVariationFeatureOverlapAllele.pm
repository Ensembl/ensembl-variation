=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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
use Bio::EnsEMBL::Variation::Utils::Constants qw(%OVERLAP_CONSEQUENCES);
use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(overlap);
use Scalar::Util qw(weaken);

our $PREDICATE_COUNT = 0;
our @SORTED_OVERLAP_CONSEQUENCES = sort {$a->tier <=> $b->tier} values %OVERLAP_CONSEQUENCES;

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
    ) = rearrange([qw(
            BASE_VARIATION_FEATURE_OVERLAP
            IS_REFERENCE
        )], @_);

    assert_ref($base_variation_feature_overlap, 'Bio::EnsEMBL::Variation::BaseVariationFeatureOverlap');
    
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
        assert_ref($bvfo, 'Bio::EnsEMBL::Variation::BaseVariationFeatureOverlap');
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

        # print STDERR "\n".join(" ", map {$_.'='.$pre->{$_}} keys %{$pre})."\n";
        
        # loop over all the possible consequences
        OC: for my $oc (@SORTED_OVERLAP_CONSEQUENCES) {
            
            # check include to hash to determine whether to run this predicate
            if(my $inc = $oc->include) {
              foreach my $key(keys %$inc) {
                if($inc->{$key} == 0) {
                  next OC if exists($pre->{$key}) && $pre->{$key} ne $inc->{$key};
                }
                else {
                  next OC unless exists($pre->{$key}) && $pre->{$key} eq $inc->{$key};
                }
              }
            }
            
            last if defined($assigned_tier) and $oc->tier > $assigned_tier;
           
            # check that this consequence applies to this type of variation feature

            if ($oc->variant_feature_class && $bvf->isa($oc->variant_feature_class)) {
                
                # check that this consequence applies to this type of feature
                if ($feat->isa($oc->feature_class)) {

                    # print STDERR $oc->SO_term." ".join(" ", map {$_.'='.$oc->include->{$_}} keys %{$oc->include})."\n";
                    
                    # if so, check if the predicate of this consequence holds for this bvfoa
                    my $check = $oc->predicate->($self, $feat, $bvfo, $bvf);
                    
                    $PREDICATE_COUNT++;
                    
                    #print STDERR $self->base_variation_feature->variation_name." ".$oc->{SO_term}." ".$self->feature->stable_id. " $check\n";

                    if ($check) {
                        push @$cons, $oc;
                        $assigned_tier = $oc->tier;
                    }
                }
            }
        }            

        $self->{overlap_consequences} = $cons;
    }
    
    return $self->{overlap_consequences};
}

# pre-calculate a bunch of attributes
# used to work out whether to execute predicates in get_all_OverlapConsequences
# uses a complex if/else structure to avoid executing unnecessary code
sub _pre_consequence_predicates {
  my ($self, $feat, $bvfo, $bvf) = @_;
  
  unless(exists($self->{pre_consequence_predicates})) {
    $bvfo ||= $self->base_variation_feature_overlap;
    $bvf  ||= $bvfo->base_variation_feature;
    $feat ||= $bvfo->feature;
    
    my $preds = {};
    
    $preds->{feature_type} = ref($feat);
    
    ## variant type
    unless(exists($bvf->{pre_consequence_predicates})) {
      
      my $bvf_preds = {};
      
      my $class_SO_term = $bvf->class_SO_term;
      my $is_sv = $bvf->isa('Bio::EnsEMBL::Variation::StructuralVariationFeature') ? 1 : 0;
      
      $bvf_preds->{sv} = $is_sv;
      
      # use SO term to determine class
      if($is_sv) {
        if($class_SO_term =~ /deletion/) {
          $bvf_preds->{deletion} = 1;
          $bvf_preds->{decrease_length} = 1;
        }
        elsif($class_SO_term =~ /insertion|duplication/) {
          $bvf_preds->{insertion} = 1;
          $bvf_preds->{increase_length} = 1;
        }
      }
      
      # otherwise for sequence variants, log the reference length here
      # we'll use it later to determine class per-allele
      else {
        $bvf_preds->{ref_length} = ($bvf->{end} - $bvf->{start}) + 1;
      }
      
      $bvf->{pre_consequence_predicates} = $bvf_preds;
    }
    
    # copy from bvf
    $preds->{$_} = $bvf->{pre_consequence_predicates}->{$_} for keys %{$bvf->{pre_consequence_predicates}};
    
    
    ## feature-specific location
    unless(exists($bvfo->{pre_consequence_predicates})) {
      
      my $bvfo_preds = {};
      my ($vf_start, $vf_end) = ($bvf->start, $bvf->end);
      my ($min_vf, $max_vf) = $vf_start > $vf_end ? ($vf_end, $vf_start) : ($vf_start, $vf_end);
      
      # within feature
      $bvfo_preds->{within_feature} = overlap($vf_start, $vf_end, $feat->start, $feat->end) ? 1 : 0;
      
      if($bvfo_preds->{within_feature}) {
        
        my ($tr_coding_start, $tr_coding_end);
        
        if($preds->{feature_type} eq 'Bio::EnsEMBL::Transcript') {
          
          # get coding region coords
          ($tr_coding_start, $tr_coding_end) = ($feat->coding_region_start, $feat->coding_region_end);
          
          # store biotype
          $bvfo_preds->{$feat->biotype} = 1;

          # use interval trees, should be whizzy fast (also neater code)
          if($bvfo->_can_use_interval_tree) {
            if(my $tree = $bvfo->_intron_interval_tree()) {
              $bvfo_preds->{intron} = scalar @{$tree->fetch($min_vf - 1, $max_vf)} ? 1 : 0;
            }

            if(my $tree = $bvfo->_intron_boundary_interval_tree()) {
              $bvfo_preds->{intron_boundary} = scalar @{$tree->fetch($min_vf - 1, $max_vf)} ? 1 : 0;
            }

            if(my $tree = $bvfo->_exon_interval_tree()) {

              # apply a "stretch" to the VF coordinates if we have a frameshift intron
              # the introns can be up to 12 bases long and if a variant falls in one
              # we actually want to call it exonic
              my $stretch = $feat->{_variation_effect_feature_cache}->{_has_frameshift_intron} ? 12 : 0;

              $bvfo_preds->{exon} = scalar @{$tree->fetch($min_vf - ($stretch + 1), $max_vf + $stretch)} ? 1 : 0;
            }
          }

          # otherwise fall back to using overlap()
          else {
            my $tr_strand = $feat->strand();
            my $vf_3_prime_end = $tr_strand > 0 ? $max_vf : $min_vf;
            
            # intronic?
            my $within_intron = 0;
            my $within_boundary = 0;
            
            foreach my $intron(@{$bvfo->_introns}) {

              my ($intron_start, $intron_end) = ($intron->start, $intron->end);

              # check intron
              if(!$within_intron && overlap($min_vf, $max_vf, $intron_start - 8, $intron_end + 8)) {
                $within_intron = 1;
              }

              if(!$within_boundary && (
                overlap($min_vf, $max_vf, $intron_start - 8, $intron_start + 8) ||
                overlap($min_vf, $max_vf, $intron_end - 8, $intron_end + 8)
              )) {
                $within_boundary = 1;
              }

              # cache whether this transcript has a frameshift intron
              # we use it later to assess exonic status
              $feat->{_variation_effect_feature_cache}->{_has_frameshift_intron} = 1 if abs($intron_end - $intron_start) <= 12;

              # we only want binary status so leave if both are true
              last if $within_intron && $within_boundary;

              # also leave if we've gone beyond the bounds of the VF
              # subsequent introns won't be in range
              last if
                ( $tr_strand > 1 && $vf_3_prime_end < ($intron->start - 8)) ||
                ( $tr_strand < 1 && $vf_3_prime_end > ($intron->end + 8));
            }

            $bvfo_preds->{intron} = $within_intron;
            $bvfo_preds->{intron_boundary} = $within_boundary;
            
            # exonic?
            $bvfo_preds->{exon} = 0;
            
            foreach my $exon(@{$bvfo->_exons}) {

              # apply a "stretch" to the exon coordinates if we have a frameshift intron
              # the introns can be up to 12 bases long and if a variant falls in one
              # we actually want to call it exonic
              my $stretch = $feat->{_variation_effect_feature_cache}->{_has_frameshift_intron} ? 12 : 0;

              # also leave if we've gone beyond the bounds of the VF
              # subsequent introns won't be in range
              last if
                ( $tr_strand > 1 && $vf_3_prime_end < $exon->start - $stretch) ||
                ( $tr_strand < 1 && $vf_3_prime_end > $exon->end + $stretch);
            
              if(overlap($min_vf, $max_vf, $exon->start - $stretch, $exon->end + $stretch)) {
                $bvfo_preds->{exon} = 1;
                last;
              }
            }
          }
        }
        
        # coding transcript specific
        if($tr_coding_start && $tr_coding_end) {

          # completely within UTR
          if(!overlap($min_vf, $max_vf, $tr_coding_start, $tr_coding_end)) {
            $bvfo_preds->{utr} = 1;
          }
        
          else {
        
            # partially within UTR
            if($min_vf < $tr_coding_start || $max_vf > $tr_coding_end) {
              $bvfo_preds->{utr} = 1;
            }

            # shortcut if trees used
            if(exists($bvfo_preds->{exon}) && !$bvfo_preds->{exon}) {
              $bvfo_preds->{non_coding} = 1;
            }

            else {
              my $cds_coords = $bvfo->cds_coords;
              
              # compeletely within coding region
              if(scalar @$cds_coords == 1 && $cds_coords->[0]->isa('Bio::EnsEMBL::Mapper::Coordinate')) {
                $bvfo_preds->{coding} = 1;
              }
              
              else {
                
                # completely in non-coding region
                if(scalar @$cds_coords == 0) {
                  $bvfo_preds->{non_coding} = 1;
                }
                
                # partly in coding region
                else {
                  $bvfo_preds->{coding} = 1;
                }
              } 
            }
          }
        }
        
        elsif($preds->{feature_type} eq 'Bio::EnsEMBL::Transcript') {
          $bvfo_preds->{non_coding} = 1;
        }
      }
      
      $bvfo->{pre_consequence_predicates} = $bvfo_preds;
    }
    
    # copy from bvfo
    $preds->{$_} = $bvfo->{pre_consequence_predicates}->{$_} for keys %{$bvfo->{pre_consequence_predicates}};
    
    
    ## allele-specific type for non-SVs
    unless($preds->{sv}) {
      my $vf_seq = $self->variation_feature_seq();
      my $ref_length = $preds->{ref_length};
    
      # most variants we see will be SNPs, so test this first
      if(length($vf_seq) == $ref_length) {
        $preds->{snp} = 1;
      }
      elsif( $ref_length > length($vf_seq) ) {
        $preds->{deletion} = 1;
        $preds->{decrease_length} = 1;
      }
      elsif( $ref_length < length($vf_seq) ) {
        $preds->{insertion} = 1;
        $preds->{increase_length} = 1;
      }
    }
    
    # copy to self
    $self->{pre_consequence_predicates} = $preds;
  }
  
  return $self->{pre_consequence_predicates};
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
    assert_ref($oc, 'Bio::EnsEMBL::Variation::OverlapConsequence');
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

