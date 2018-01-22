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

Bio::EnsEMBL::Variation::BaseVariationFeatureOverlap

=head1 SYNOPSIS

    use Bio::EnsEMBL::Variation::BaseVariationFeatureOverlap;

    my $bvfo = Bio::EnsEMBL::Variation::BaseVariationFeatureOverlap->new(
        -feature                => $feature,
        -base_variation_feature => $var_feat
    );

    print "consequence type: ", (join ",", @{ $bvfo->consequence_type }), "\n";
    print "most severe consequence: ", $bvfo->display_consequence, "\n";

=head1 DESCRIPTION

A BaseVariationFeatureOverlap represents a BaseVariationFeature which is in close
proximity to another Ensembl Feature. It is the superclass of variation feature
specific classes such as VariationFeatureOverlap and StructuralVariationOverlap
and has methods common to all such objects. You will not normally instantiate this
class directly, instead instantiating one of the more specific subclasses.

=cut

package Bio::EnsEMBL::Variation::BaseVariationFeatureOverlap;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(overlap within_cds);

sub new {
  my $class = shift;

  my (
    $adaptor,
    $base_variation_feature,
    $feature, 
    $no_transfer
  );

  if($Bio::EnsEMBL::Utils::Argument::NO_REARRANGE) {
    my %args = @_;

    (
      $adaptor,
      $base_variation_feature,
      $feature, 
      $no_transfer
    ) = (
      $args{-adaptor},
      $args{-base_variation_feature},
      $args{-feature},
      $args{-no_transfer}
    );
  }

  else {
    (
      $adaptor,
      $base_variation_feature,
      $feature, 
      $no_transfer
    ) = rearrange([qw(
      ADAPTOR
      BASE_VARIATION_FEATURE
      FEATURE
      NO_TRANSFER
    )], @_);
  }

  assert_ref($base_variation_feature, 'Bio::EnsEMBL::Variation::BaseVariationFeature') if $Bio::EnsEMBL::Utils::Scalar::ASSERTIONS;

  if (defined $feature) {

    assert_ref($feature, 'Bio::EnsEMBL::Feature') if $Bio::EnsEMBL::Utils::Scalar::ASSERTIONS;

    # we need to ensure the Feature and the BaseVariationFeature live on the same slice
    # so we explicitly transfer the Feature here
    unless($no_transfer && $no_transfer == 1) {
      $feature = $feature->transfer($base_variation_feature->slice) 
        or throw("Unable to transfer the supplied feature to the same slice as the base variation feature");
    }
  }

  my $self = bless {
    base_variation_feature  => $base_variation_feature,
    feature                 => $feature,
    adaptor                 => $adaptor,
  }, $class;

  return $self;
}

sub new_fast {
    my ($class, $hashref) = @_;
    return bless $hashref, $class;
}

=head2 feature

  Arg [1]    : (optional) A Bio::EnsEMBL::Feature
  Description: Get/set the associated Feature, lazy-loading it if required
  Returntype : Bio::EnsEMBL::Feature
  Exceptions : throws isf the argument is the wrong type
  Status     : Stable

=cut

sub feature {
    my ($self, $feature, $type) = @_;
    
    if ($feature) {
        assert_ref($feature, 'Bio::EnsEMBL::Feature') if $Bio::EnsEMBL::Utils::Scalar::ASSERTIONS;
        $self->{feature} = $feature;
    }
 
    if ($type && !$self->{feature}) {
        if (my $adap = $self->{adaptor}) {
            # try to lazy load the feature
            my $feature_stable_id = $self->{_feature_stable_id};
            if ($type eq 'RegulatoryFeature') {
                my $funcgen_db = $adap->db->get_db_adaptor('funcgen');
                unless ($funcgen_db) {
                    warn("Ensembl Funcgen DB is missing.");
                }
                my $rfa = $funcgen_db->get_RegulatoryFeatureAdaptor;
                if ($feature_stable_id) {
                    my $feature = $rfa->fetch_by_stable_id($feature_stable_id);
                    $self->{feature} = $feature;
                }
            } elsif ($type eq 'MotifFeature') {
                my $funcgen_db = $adap->db->get_db_adaptor('funcgen');
                unless ($funcgen_db) {
                    warn("Ensembl Funcgen DB is missing.")
                }
                my $mfa = $funcgen_db->get_MotifFeatureAdaptor;
                if (my $motif_feature_id = $self->{motif_feature_id}) {
                    my $feature = $mfa->fetch_by_dbID($motif_feature_id);
                    $self->{feature} = $feature;
                }    
            } else {
                my $get_method = 'get_'.$type.'Adaptor';
               
                # XXX: this can doesn't work because the method is AUTOLOADed, need to rethink this...
                #if ($adap->db->dnadb->can($get_method)) {
                    if (my $fa = $adap->db->dnadb->$get_method) {
                        
                        # if we have a stable id for the feature use that
                        if (my $feature_stable_id = $self->{_feature_stable_id}) {
                            if (my $f = $fa->fetch_by_stable_id($feature_stable_id)) {
                                $self->{feature} = $f;
                                #delete $self->{_feature_stable_id};
                            }
                        }
                        elsif (my $feature_label = $self->{_feature_label}) {
                            # get a slice covering the vf
                            
                            #for my $f ($fa->fetch_all_by_Slice_constraint)
                        }
                    }
                #}
                else {
                    warn "Cannot get an adaptor for type: $type";
                }
            }
        }
    }
    
    return $self->{feature};
}

sub _fetch_feature_for_stable_id {
    
    # we shouldn't actually need this method as there will apparently
    # soon be core support for fetching any feature by its stable id, 
    # but I'm waiting for core to add this...

    my ($self, $feature_stable_id) = @_;
    
    my $type_lookup = {
        G   => { type => 'Gene',                 group => 'core' },
        T   => { type => 'Transcript',           group => 'core'  },
        R   => { type => 'RegulatoryFeature',    group => 'funcgen' },
    };
    
    if ($feature_stable_id =~ /^ENS[A-Z]*([G|R|T])\d+$/) {
        
        my $type  = $type_lookup->{$1}->{type};
        my $group = $type_lookup->{$1}->{group};
        
        if (my $adap = $self->{adaptor}) {
            
            my $get_method = 'get_'.$type.'Adaptor';
            
            if ($adap->db->dnadb->can($get_method)) {
                if (my $fa = $adap->db->dnadb->$get_method) {
                    
                    # if we have a stable id for the feature use that
                    if (my $feature_stable_id = $self->{_feature_stable_id}) {
                        if (my $f = $fa->fetch_by_stable_id($feature_stable_id)) {
                            $self->{feature} = $f;
                            #delete $self->{_feature_stable_id};
                        }
                    }
                    elsif (my $feature_label = $self->{_feature_label}) {
                        # get a slice covering the vf
                        
                        
                        #for my $f ($fa->fetch_all_by_Slice_constraint)
                    }
                }
            }
            else {
                warn "Cannot get an adaptor for type: $type";
            }
    }
    }
}

sub _fetch_adaptor_for_group {
    my ($self, $group) = @_;
    
}

sub _feature_stable_id {
    my $self = shift;
    
    if (my $id = $self->{_feature_stable_id}) {
        return $id;
    }
    elsif ($self->feature && $self->feature->can('stable_id')) {
        $self->{_feature_stable_id} = $self->feature->stable_id;
        
        return $self->{_feature_stable_id};
    }
    else {
        return undef;
    }
}

=head2 base_variation_feature

  Arg [1]    : (optional) A Bio::EnsEMBL::Variation::BaseVariationFeature
  Description: Get/set the associated BaseVariationFeature
  Returntype : Bio::EnsEMBL::Variation::BaseVariationFeature
  Exceptions : throws if the argument is the wrong type
  Status     : Stable

=cut

sub base_variation_feature {
    my ($self, $bvf) = @_;
    
    if ($bvf) {
        assert_ref($bvf, 'Bio::EnsEMBL::Variation::BaseVariationFeature') if $Bio::EnsEMBL::Utils::Scalar::ASSERTIONS;
        $self->{base_variation_feature} = $bvf;
    }

    return $self->{base_variation_feature};
}

=head2 add_BaseVariationFeatureOverlapAllele

  Arg [1]    : A Bio::EnsEMBL::Variation::BaseVariationFeatureOverlapAllele instance
  Description: Add an allele to this BaseVariationFeatureOverlap
  Returntype : none
  Exceptions : throws if the argument is not the expected type
  Status     : At Risk

=cut

sub add_BaseVariationFeatureOverlapAllele {
    my ($self, $bvfoa) = @_;

    assert_ref($bvfoa, 'Bio::EnsEMBL::Variation::BaseVariationFeatureOverlapAllele') if $Bio::EnsEMBL::Utils::Scalar::ASSERTIONS;

    if ($bvfoa->is_reference) {
        $self->{reference_allele} = $bvfoa;
    }
    else {
        my $alt_alleles = $self->{alt_alleles} ||= [];
        push @$alt_alleles, $bvfoa;
    }
}

=head2 get_reference_BaseVariationFeatureOverlapAllele

  Description: Get the object representing the reference allele of this BaseVariationFeatureOverlapAllele
  Returntype : Bio::EnsEMBL::Variation::BaseVariationFeatureOverlapAllele instance
  Exceptions : none
  Status     : At Risk

=cut

sub get_reference_BaseVariationFeatureOverlapAllele {
    my $self = shift;
    return $self->{reference_allele};
}

=head2 get_all_alternate_BaseVariationFeatureOverlapAlleles

  Description: Get a list of the alternate alleles of this BaseVariationFeatureOverlapAllele
  Returntype : listref of Bio::EnsEMBL::Variation::BaseVariationFeatureOverlapAllele objects
  Exceptions : none
  Status     : At Risk

=cut

sub get_all_alternate_BaseVariationFeatureOverlapAlleles {
    my $self = shift;

    $self->{alt_alleles} ||= [];
    
    return $self->{alt_alleles};
}

=head2 get_all_BaseVariationFeatureOverlapAlleles

  Description: Get a list of the all the alleles, both reference and alternate, of this
               BaseVariationFeatureOverlap
  Returntype : listref of Bio::EnsEMBL::Variation::BaseVariationFeatureOverlapAllele objects
  Exceptions : none
  Status     : At Risk

=cut

sub get_all_BaseVariationFeatureOverlapAlleles {
    my $self = shift;
	
    my @alleles = @{ $self->get_all_alternate_BaseVariationFeatureOverlapAlleles };
	
    my $ref = $self->get_reference_BaseVariationFeatureOverlapAllele;
	
    unshift @alleles, $ref if defined $ref;
	
    return \@alleles;
}

=head2 consequence_type

  Arg [1]    : (optional) String $term_type
  Description: Get a list of all the unique consequence terms of the alleles of this 
               BaseVariationFeatureOverlap. By default returns Ensembl display terms
               (e.g. 'NON_SYNONYMOUS_CODING'). $term_type can also be 'label'
               (e.g. 'Non-synonymous coding'), 'SO' (Sequence Ontology, e.g.
               'non_synonymous_codon') or 'NCBI' (e.g. 'missense')
  Returntype : listref of strings
  Exceptions : none
  Status     : Stable

=cut

sub consequence_type {
    my $self = shift;
	my $term_type = shift;
	
	my $method_name;
	
    # delete cached term
    if(defined($term_type)) {
        delete $self->{_consequence_type};
		$method_name = $term_type.($term_type eq 'label' ? '' : '_term');
		$method_name = 'SO_term' unless defined $self->most_severe_OverlapConsequence && $self->most_severe_OverlapConsequence->can($method_name);
    }
	
	$method_name ||= 'SO_term';
    
    unless ($self->{_consequence_type}) {
        
        # use a hash to ensure we don't include redundant terms (because more than one
        # allele may have the same consequence SO_term)

        my %cons_types;

        for my $allele (@{ $self->get_all_alternate_BaseVariationFeatureOverlapAlleles }) {
            for my $cons (@{ $allele->get_all_OverlapConsequences }) {
                $cons_types{$cons->$method_name} = $cons->rank;
            }
        }

        # sort the consequence types by rank such that the more severe terms are earlier in the list

        $self->{_consequence_type} = [ sort { $cons_types{$a} <=> $cons_types{$b} } keys %cons_types ];
    }
    
    return $self->{_consequence_type};
}

=head2 most_severe_OverlapConsequence

  Description: Get the OverlapConsequence considered (by Ensembl) to be the most severe 
               consequence of all the alleles of this VariationFeatureOverlap 
  Returntype : Bio::EnsEMBL::Variation::OverlapConsequence
  Exceptions : none
  Status     : At Risk

=cut

sub most_severe_OverlapConsequence {
    my $self = shift;
    
    unless ($self->{_most_severe_consequence}) {
        
        my $highest;
        
        for my $allele (@{ $self->get_all_alternate_BaseVariationFeatureOverlapAlleles }) {
            for my $cons (@{ $allele->get_all_OverlapConsequences }) {
                $highest ||= $cons;
                if ($cons->rank < $highest->rank) {
                    $highest = $cons;
                }
            }
        }
        
        $self->{_most_severe_consequence} = $highest;
    }
    
    return $self->{_most_severe_consequence};
}

=head2 display_consequence

  Arg [1]    : (optional) String $term_type
  Description: Get the term for the most severe OverlapConsequence of this 
               VariationFeatureOverlap. By default returns Ensembl display terms
               (e.g. 'NON_SYNONYMOUS_CODING'). $term_type can also be 'label'
               (e.g. 'Non-synonymous coding'), 'SO' (Sequence Ontology, e.g.
               'non_synonymous_codon') or 'NCBI' (e.g. 'missense')
  Returntype : string
  Exceptions : none
  Status     : Stable

=cut

sub display_consequence {
    my $self = shift;
	my $term_type = shift;
	
	my $method_name;
	
    # delete cached term
    if(defined($term_type)) {
		$method_name = $term_type.($term_type eq 'label' ? '' : '_term');
		$method_name = 'SO_term' unless @{$self->get_all_OverlapConsequences} && $self->get_all_OverlapConsequences->[0]->can($method_name);
    }
	
	$method_name ||= 'SO_term';
    
    my $worst_conseq = $self->most_severe_OverlapConsequence;

    return $worst_conseq ? $worst_conseq->$method_name : '';
}

sub adaptor {
    my $self = shift;
    $self->{adaptor} = shift if @_;
    
    # make adaptor an anonymous hash in its absence
    # this allows the VEP to cache OverlapConsequences in offline mode
    $self->{adaptor} ||= {};
    
    return $self->{adaptor};
}

1;

