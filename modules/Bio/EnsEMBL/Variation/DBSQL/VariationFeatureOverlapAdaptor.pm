=head1 LICENSE

 Copyright (c) 1999-2012 The European Bioinformatics Institute and
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

Bio::EnsEMBL::Variation::DBSQL::VariationFeatureOverlapAdaptor

=head1 DESCRIPTION

This is the superclass of all Adaptors that fetch VariationFeatureOverlap
objects and their various subclasses, and it provides methods common to
all such adaptors, such as fetching by VariationFeature. You should not
generally use this class directly, but instead use one of the feature
specific adaptors such as the TranscriptVariationAdaptor.

=cut
 
use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::VariationFeatureOverlapAdaptor;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use base qw(Bio::EnsEMBL::Variation::DBSQL::BaseAdaptor);

sub new_fake {
    my $class = shift;
    my $species = shift;

    my $self = bless {}, $class;

    return $self;
}

=head2 fetch_all_by_Features

  Arg [1]    : listref of Bio::EnsEMBL::Features, or subclasses
  Description: Fetch all germline VariationFeatureOverlap objects associated 
               with the given list of Features
  Returntype : listref of Bio::EnsEMBL::Variation::VariationFeatureOverlap objects
  Status     : At risk

=cut

sub fetch_all_by_Features {
    my ($self, $features) = @_;
    return $self->fetch_all_by_Features_with_constraint($features,'somatic = 0');
}

=head2 fetch_all_somatic_by_Features

  Arg [1]    : listref of Bio::EnsEMBL::Features, or subclasses
  Description: Fetch all somatic VariationFeatureOverlap objects associated 
               with the given list of Features
  Returntype : listref of Bio::EnsEMBL::Variation::VariationFeatureOverlap objects
  Status     : At risk

=cut

sub fetch_all_somatic_by_Features {
    my ($self, $features) = @_;
    return $self->fetch_all_by_Features_with_constraint($features,'somatic = 1');
}

=head2 fetch_all_by_Features_with_constraint

  Arg [1]    : listref of Bio::EnsEMBL::Features, or subclasses
  Arg [2]    : extra SQL constraint for the query 
  Description: Fetch all VariationFeatureOverlap objects associated 
               with the given list of Features
  Returntype : listref of Bio::EnsEMBL::Variation::VariationFeatureOverlap objects
  Status     : At risk

=cut

sub fetch_all_by_Features_with_constraint {
    
    my ($self, $features, $constraint) = @_;
   
    my %feats_by_id = map { $_->stable_id => $_ } @$features;
    
    my $id_str = join ',', map {"'$_'"} keys %feats_by_id;
    
    my $full_constraint = "feature_stable_id in ( $id_str )";
    $full_constraint .= " AND $constraint" if $constraint;
    
    my $vfos = $self->generic_fetch($full_constraint);
    
    for my $vfo (@$vfos) {
        if ($vfo->{_feature_stable_id}) {
            my $feat_id = delete $vfo->{_feature_stable_id};
            $vfo->{feature} = $feats_by_id{$feat_id};
        }
    }
    
    return $vfos;
}

=head2 fetch_all_by_VariationFeatures

  Arg [1]    : listref of Bio::EnsEMBL::Variation::VariationFeatures
  Arg [2]    : (optional) listref of Bio::EnsEMBL::Features to further limit the query
  Description: Fetch all VariationFeatureOverlap objects associated 
               with the given list of VariationFeatures
  Returntype : listref of Bio::EnsEMBL::Variation::VariationFeatureOverlap objects
  Status     : At risk

=cut

sub fetch_all_by_VariationFeatures {
    
    my ($self, $vfs, $features) = @_;
   
    my %vfs_by_id = map { $_->dbID => $_ } @$vfs;
  
    my $id_str = join ',', keys %vfs_by_id;

    my $constraint = "variation_feature_id in ( $id_str )";

    my $vfos;

    if ($features) {
        # if we're passed some features, fetch by features with the VF ids as an 
        # extra constraint
        $vfos = $self->fetch_all_by_Features_with_constraint($features, $constraint);
    }
    else {
        # otherwise just fetch the VFs directly
        $vfos = $self->generic_fetch($constraint);
    }

    # attach the VariationFeatures to the VariationFeatureOverlaps because we have them already

    for my $vfo (@$vfos) {
        if ($vfo->{_variation_feature_id}) {
            $vfo->variation_feature($vfs_by_id{delete $vfo->{_variation_feature_id}});
        }
    }
   
    return $vfos;
}

sub _get_term_object {
    my ($self, $term) = @_;

    my $oa = $self->{_ontology_adaptor} ||= 
        Bio::EnsEMBL::Registry->get_adaptor( 'Multi', 'Ontology', 'OntologyTerm' );

    my $terms = $oa->fetch_all_by_name($term, 'SO');

    if (@$terms > 1) {
        warn "Ambiguous term '$term', just using first result";
    }

    return $terms->[0];
}

sub _get_child_terms {
    my ($self, $parent_term) = @_;

    my $parent_obj = $self->_get_term_object($parent_term);

    my $all_terms = $parent_obj->descendants;

    unshift @$all_terms, $parent_obj;

    return $all_terms;
}

sub _get_parent_terms {
    my ($self, $child_term) = @_;

    my $child_obj = $self->_get_term_object($child_term);

    my $all_terms = $child_obj->ancestors;

    unshift @$all_terms, $child_obj;

    return $all_terms;
}

sub _get_VariationFeatureOverlapAlleles_under_SO_term {
    my ($self, $term, $vfoas) = @_;

    my $terms = $self->_get_child_terms($term);

    my @found;

    ALLELES : for my $vfoa (@$vfoas) {
        for my $cons (@{ $vfoa->get_all_OverlapConsequences }) {
            for my $term (@$terms) {
                if ($cons->SO_term eq $term->name) {
                    push @found, $vfoa;
                    next ALLELES;
                }
            }
        }
    }

    return \@found;
}

sub _get_consequence_constraint {
    
    my ($self, $query_term) = @_;

    # we allow either an ontology term object, or just a string
    $query_term = UNIVERSAL::can($query_term, 'name') ? $query_term->name : $query_term;
   
    # get a hash mapping consequence terms to numerical values (specifically powers of 2)
    my $cons_map = $self->_consequence_type_map('transcript_variation', 'consequence_types');

    # we store only the most specific consequence term, so we need to get all children of 
    # the query term
    my $terms = $self->_get_child_terms($query_term);

    # and then build up the numerical value for our query by ORing together all the terms
    my $query = 0;

    for my $term (@$terms) {
        next unless $cons_map->{$term->name};
        $query |= $cons_map->{$term->name};
    }

    unless ($self->{_possible_consequences}) {

        # we need a list of the numerical values of all possible 
        # consequence term combinations we have actually observed

        my $sth = $self->dbc->prepare(qq{
            SELECT DISTINCT(consequence_types)
            FROM transcript_variation
        });

        $sth->execute;

        my $cons;

        $sth->bind_columns(\$cons);

        my @poss_cons;

        while ($sth->fetch) {
            # construct the numerical value by ORing together each combination
            # (this is much quicker than SELECTing consequence_types+0 above which
            # is what I used to do, but this seems to mean the db can't use the index)
        
            my $bit_val = 0;
            
            for my $term (split /,/, $cons) {
                $bit_val |= $cons_map->{$term};
            }

            push @poss_cons, $bit_val;
        }

        $self->{_possible_consequences} = \@poss_cons;
    }

    # we now find all combinations that include our query by ANDing 
    # the query with all possible combinations and combine these into 
    # our query string

    my $id_str = join ',', grep { $_ & $query } @{ $self->{_possible_consequences} }; 

    my $constraint = "consequence_types IN ($id_str)"; 

    return $constraint;
}

sub fetch_all_by_SO_term {
    my ($self, $term) = @_;

    my $constraint = $self->_get_consequence_constraint($term);

    return $self->generic_fetch($constraint);
}


1;
