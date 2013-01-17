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
  Status     : Stable

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
  Status     : Stable

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
  Status     : Stable

=cut

sub fetch_all_by_Features_with_constraint {
    my $self = shift;
    
    my ($features, $constraint) = @_;
   
    my $vfos = $self->_func_all_by_Features_with_constraint(@_, 'fetch');
    
    # Note duplicated code 
    my %feats_by_id = map { $_->stable_id => $_ } @$features;

    for my $vfo (@$vfos) {
        if ($vfo->{_feature_stable_id}) {
            my $feat_id = delete $vfo->{_feature_stable_id};
            $vfo->{feature} = $feats_by_id{$feat_id};
        }
    }
    
    return $vfos;
}

sub _func_all_by_Features_with_constraint {
    my ($self, $features, $constraint, $func) = @_;
   
    my %feats_by_id = map { $_->stable_id => $_ } @$features;
    
    my $id_str = join ',', map {"'$_'"} keys %feats_by_id;
    
    my $full_constraint = "feature_stable_id in ( $id_str )";
    $full_constraint .= " AND $constraint" if $constraint;

    my $method = "generic_" . $func;
    my $data = $self->$method($full_constraint);

    return $data;
}

sub count_all_by_Features_with_constraint {
    my $self = shift;
    my ($features, $constraint) = @_;

    my $count = $self->_func_all_by_Features_with_constraint(@_, 'count');

    if (!defined($count)) { $count = 0; }
 
    return $count;
}

=head2 fetch_all_by_VariationFeatures

  Arg [1]    : listref of Bio::EnsEMBL::Variation::VariationFeatures
  Arg [2]    : (optional) listref of Bio::EnsEMBL::Features to further limit the query
  Description: Fetch all VariationFeatureOverlap objects associated 
               with the given list of VariationFeatures
  Returntype : listref of Bio::EnsEMBL::Variation::VariationFeatureOverlap objects
  Status     : Stable

=cut
sub fetch_all_by_VariationFeatures {
    my ($self, $vfs, $features) = @_;
    return $self->fetch_all_by_VariationFeatures_with_constraint($vfs,$features,undef);
}
    
sub count_all_by_VariationFeatures {
    my ($self, $vfs, $features) = @_;
    return $self->count_all_by_VariationFeatures_with_constraint($vfs,$features,undef);
}

sub count_all_by_VariationFeatures_with_constraint {
    my $self = shift;
    my ($vfs, $features, $constraint) = @_;

    my $allcounts = $self->_func_all_by_VariationFeatures_with_constraint(@_ , 'count');

    my $total = 0;
    for my $count (@$allcounts) {
        $total += $count;
    }
   
    return $total;
}

sub _func_all_by_VariationFeatures_with_constraint {
    my ($self, $vfs, $features, $constraint, $func) = @_;

    my %vfs_by_id = map { $_->dbID => $_ } @$vfs;

    my @vfids = keys %vfs_by_id;

    if (!scalar(@vfids)) {
      return [];
    }

    my @alldata;

    while (@vfids) {
  
      my $fullconstraint = $constraint;

      my @vfid_subset = splice(@vfids,0,50000);

      my $id_str = join ',', @vfid_subset;
  
      if ($id_str eq '') {
        last;
      }
  
      if ($fullconstraint) {
        $fullconstraint .= " AND ";
      }
      $fullconstraint .= "variation_feature_id in ( $id_str )";
  
  
      my $data;
  
      if ($features) {
          # if we're passed some features, fetch/count by features with the VF ids as an 
          # extra constraint
          my $method = $func . "_all_by_Features_with_constraint";
          $data = $self->$method($features, $fullconstraint);
      }
      else {
          # otherwise just fetch/count the VFs directly
          my $method = "generic_" . $func;
          $data = $self->$method($fullconstraint);
      }
      push @alldata,ref($data) eq 'ARRAY' ? @$data : $data;
    }

    return \@alldata;
} 

sub fetch_all_by_VariationFeatures_with_constraint {
    my $self = shift;
    my ($vfs, $features, $constraint) = @_;

    my $allvfos = $self->_func_all_by_VariationFeatures_with_constraint(@_ , 'fetch');
   

    my %vfs_by_id = map { $_->dbID => $_ } @$vfs;

    # attach the VariationFeatures to the VariationFeatureOverlaps because we have them already

    for my $vfo (@$allvfos) {
        if ($vfo->{_variation_feature_id}) {
            $vfo->variation_feature($vfs_by_id{delete $vfo->{_variation_feature_id}});
        }
    }
   
    return $allvfos;
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

# call to method in BaseAdaptor
sub _get_consequence_constraint {
    my $self = shift;
	return $self->SUPER::_get_consequence_constraint('transcript_variation', @_);
}

sub fetch_all_by_SO_terms {
    my ($self, $terms) = @_;

    my $constraint = $self->_get_consequence_constraint($terms);

    return $self->generic_fetch($constraint);
}


1;
