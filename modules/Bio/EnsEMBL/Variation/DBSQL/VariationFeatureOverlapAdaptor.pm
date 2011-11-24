=head1 LICENSE

 Copyright (c) 1999-2011 The European Bioinformatics Institute and
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

1;
