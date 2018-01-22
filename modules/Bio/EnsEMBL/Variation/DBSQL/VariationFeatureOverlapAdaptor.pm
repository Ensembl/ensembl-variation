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
use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(MAX_DISTANCE_FROM_TRANSCRIPT);
use Scalar::Util qw(looks_like_number);

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

    my $use_vcf = $self->db->use_vcf();
    my @vcf_vfos;

    if($use_vcf) {
      my $vca = $self->db->get_VCFCollectionAdaptor;
      my @vcfs = grep {$_->use_as_source} @{$vca->fetch_all};

      foreach my $f(@$features) {
        my $expanded_slice = $f->feature_Slice->expand(MAX_DISTANCE_FROM_TRANSCRIPT, MAX_DISTANCE_FROM_TRANSCRIPT);
        my $f_slice        = $f->slice;

        my @vfs =
          map {$_->transfer($f_slice)}
          map {@{$_->get_all_VariationFeatures_by_Slice($expanded_slice, 1)}}
          @vcfs;

        my $strong_ref_copy;
        push @vcf_vfos,
          map {$strong_ref_copy = $_->{base_variation_feature}; $_->{base_variation_feature} = $strong_ref_copy; $_}
          @{$self->_fetch_all_by_VariationFeatures_no_DB(\@vfs, [$f], undef, 1)};
      }

      if($use_vcf == 2) {
        return $func eq 'count' ? scalar @vcf_vfos : \@vcf_vfos;
      }
    }
   
    my %feats_by_id = map { $_->stable_id => $_ } @$features;
    
    my $id_str = join ',', map {"'$_'"} keys %feats_by_id;
    
    my $full_constraint = "feature_stable_id in ( $id_str )";
    $full_constraint .= " AND $constraint" if $constraint;

    my $method = "generic_" . $func;
    my $data = $self->$method($full_constraint);

    if($func eq 'count') {
      return (scalar @vcf_vfos + $data);
    }
    else {
      push @vcf_vfos, @$data;
      return \@vcf_vfos;
    }
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
    
    # split into those with a real dbID and those without
    my (@with_id, @no_id);
    foreach my $vf(@$vfs) {
      if(looks_like_number($vf->dbID)) {
        push @with_id, $vf;
      }
      else {
        push @no_id, $vf;
      }
    }

    my @alldata;
    
    # deal with those with no ID
    if(scalar @no_id) {
      my $method = '_'.$func.'_all_by_VariationFeatures_no_DB';

      # try to fetch VFOs from VFs
      # probably not reliable to do this unless we have $features
      # it's only really going to be of benefit for Transcript lookups too
      my $remaining = [];

      if($features && @$features && ref($features->[0]) eq 'Bio::EnsEMBL::Transcript') {
        while(my $vf = shift @no_id) {
          if(my $vfo_hash = $vf->{transcript_variations}) {
            my @tmp_vfos = grep {defined($_)} map {$vfo_hash->{$vf->_get_transcript_key($_)}} @$features;

            if(scalar @tmp_vfos == scalar @$features) {
              push @alldata, $func eq 'count' ? scalar @tmp_vfos : @tmp_vfos;
            }
            else {
              push @$remaining, $vf;
            }
          }
          else {
            push @$remaining, $vf;
          }          
        }
      }
      else {
        $remaining = \@no_id;
      }

      if(@$remaining) {
        $_->reset_consequence_data for @$remaining;
        my $data = $self->$method($remaining, $features, $constraint);
        push @alldata, ref($data) eq 'ARRAY' ? @$data : $data;
      }
    }

    my %vfs_by_id = map { $_->dbID => $_ } grep {$_->dbID} @with_id;

    my @vfids = keys %vfs_by_id;

    if(!scalar(@vfids)) {
      return \@alldata;
    }

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
   

    my %vfs_by_id = map { $_->dbID => $_ } grep {$_->dbID} @$vfs;

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

## this method fetches ranged slices from variation features
## the idea is to create as small a number of slices as possible from
## an arrayref of VFs
sub _get_ranged_slices_from_VariationFeatures {
  my ($self, $vfs, $range_size) = @_;

  return [] unless @$vfs;

  # quick check - do all the VFs have the same slice?
  my %slice_refs = map {$_->slice + 0 => 1} @$vfs;
  if(scalar keys %slice_refs == 1) {
    return [$vfs->[0]->slice->expand(MAX_DISTANCE_FROM_TRANSCRIPT, MAX_DISTANCE_FROM_TRANSCRIPT)];
  }

  # default range
  $range_size ||= MAX_DISTANCE_FROM_TRANSCRIPT;

  # create a set of non-overlapping ranges covering all variants
  my @ranges;
  my ($prev_start, $prev_chr);

  foreach my $vf(
    sort {
      $_->seq_region_name cmp $_->seq_region_name ||
      $_->seq_region_start <=> $_->seq_region_end
    } @$vfs
  ) {
    my $vf_start = $vf->seq_region_start;
    my $vf_chr   = $vf->seq_region_name;

    # start a new range if new chromosome or variant beyond previous range
    if(
      !$prev_start || !$prev_chr ||
      $prev_chr ne $vf_chr ||
      $vf_start - $prev_start > $range_size
    ) {
      push @ranges, {
        id => $vf->slice->get_seq_region_id,
        s  => $vf_start - $range_size,
        e  => $vf_start + $range_size,
      };
    }

    # expand previous range if var falls within it
    elsif(@ranges) {
      $ranges[-1]->{e} = $vf_start + $range_size;
    }

    $prev_chr   = $vf_chr;
    $prev_start = $vf_start;
  }

  # convert ranges to slices and return
  my $sa = $vfs->[0]->slice->adaptor;

  return [
    map {
      $sa->fetch_by_seq_region_id(
        $_->{id}, $_->{s}, $_->{e}
      )
    } @ranges
  ];
}


1;
