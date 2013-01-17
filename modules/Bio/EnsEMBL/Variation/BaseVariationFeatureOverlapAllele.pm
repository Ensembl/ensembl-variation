=head1 LICENSE

 Copyright (c) 1999-2013 The European Bioinformatics Institute and
 Genome Research Limited.  All rights reserved.

 This software is distributed under a modified Apache license.
 For license details, please see

   http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <dev@ensembl.org>.

 Questions may also be sent to the Ensembl help desk at
 <helpdesk@ensembl.org>.

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
use Scalar::Util qw(weaken);

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
    my ($class, $hashref) = @_;
    my $self = bless $hashref, $class;
    # avoid a memory leak, because the bvfo also has a reference to us
    weaken $self->{base_variation_feature_overlap} if $self->{base_variation_feature_overlap};
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
        
        # loop over all the possible consequences
        for my $oc (@{$self->get_sorted_OverlapConsequences}) {
            
            last if defined($assigned_tier) and $oc->tier > $assigned_tier;
           
            # check that this consequence applies to this type of variation feature

            if ($oc->variant_feature_class && $self->base_variation_feature->isa($oc->variant_feature_class)) {
                
                # check that this consequence applies to this type of feature

                if ($self->feature->isa($oc->feature_class)) {
                    
                    # if so, check if the predicate of this consequence holds for this bvfoa
                    my $check = $oc->predicate->($self);
                    
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

sub get_sorted_OverlapConsequences {
    my $self = shift;
    
    if(!defined($self->base_variation_feature_overlap->adaptor->{sorted_cons})) {
        my @sorted = sort {$a->tier <=> $b->tier} values %OVERLAP_CONSEQUENCES;
        $self->base_variation_feature_overlap->adaptor->{sorted_cons} = \@sorted;
    }
    
    return $self->base_variation_feature_overlap->adaptor->{sorted_cons};
}

1;

