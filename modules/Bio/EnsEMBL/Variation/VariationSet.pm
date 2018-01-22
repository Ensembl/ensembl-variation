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

# Ensembl module for Bio::EnsEMBL::Variation::VariationSet
#
#

=head1 NAME

Bio::EnsEMBL::Variation::VariationSet - Ensembl representation of a set of
variations.

=head1 SYNOPSIS

  use Bio::EnsEMBL::Variation::VariationSet;

  ...



=head1 DESCRIPTION

This is a class representing a set of variations that are grouped by e.g.
study, method, quality measure etc.

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::Variation::VariationSet;

use Bio::EnsEMBL::Storable;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Iterator;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Storable);


=head2 new

  Arg [dbID] :
    int - unique internal identifier for this allele group

  Arg [ADAPTOR] :
    Bio::EnsEMBL::Variation::DBSQL::VariationSetAdaptor

  Arg [NAME] :
    string - the name of this variation set

  Arg [DESCRIPTION] :
    string - A description explaining the charcteristics of this variation set

  Arg [SHORT_NAME] :
    string - the short name of this variation set

  Example    :
    $ag = Bio::EnsEMBL::Variation::VariationSet->new
      (
       -dbID => 12,
       -adaptor => $var_set_adaptor,
       -name   => 'Phenotype-associated variations',
       -description => 'Variations that have been associated with a phenotype',
       -short_name => 'ph_variants'
      );
  Description: Constructor.  Instantiates a new VariationSet
  Returntype : Bio::EnsEMBL::Variation::VariationSet
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub new {
  my $class = shift;

  my ($dbID, $adaptor, $name, $description, $short_name) =
    rearrange([qw(DBID ADAPTOR NAME DESCRIPTION SHORT_NAME)], @_);
  
  # Check that the dbID does not exceed the maximum dbID that can be stored in the variation_set_id SET construct in variation_set_variation
  warn("Primary key variation_set_id $dbID for variation set '$name' exceeds " . $Bio::EnsEMBL::Variation::DBSQL::VariationSetAdaptor::MAX_VARIATION_SET_ID . ". Therefore, this variation set cannot be properly referenced in variation_set_variation") if ($dbID > $Bio::EnsEMBL::Variation::DBSQL::VariationSetAdaptor::MAX_VARIATION_SET_ID);
  
  return bless {'dbID' => $dbID,
                'adaptor' => $adaptor,
                'name' => $name,
                'description' => $description,
                'short_name' => $short_name}, $class;
}

=head2 description

  Arg [1]    : string $description
  Example    : print $vs->description();
  Description: Getter/Setter for the description of this VariationSet
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub description {
  my $self = shift;
  my $desc = shift;
  
  $self->{'description'} = $desc if (defined($desc));
  
  return $self->{'description'};
}

=head2 get_all_sub_VariationSets
  Arg [1]    : (optional) boolean $only_immediate
               If true, will only get the direct subsets of this variation. The default behaviour is
               to recursively get all subsets.
  Example    : print $vs->get_all_sub_VariationSets();
  Description: Recursively gets all variation sets that are subsets of this variation set.
  Returntype : reference to list of Bio::EnsEMBL::Variation::VariationSet
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_sub_VariationSets {
  my $self = shift;
  my $only_immediate = shift;
  
  # A database adaptor must be attached to this object   
  if(!$self->adaptor()) {
    warning('Cannot get sub variation sets without attached adaptor');
    return [];
  }
  
  return $self->adaptor->fetch_all_by_super_VariationSet($self,$only_immediate);
}

=head2 get_all_super_VariationSets
  Arg [1]    : (optional) boolean $only_immediate
               If true, will only get the direct supersets of this variation. The default behaviour is
               to recursively get all supersets.
  Example    : print $vs->get_all_super_VariationSets();
  Description: Recursively gets all variation sets that are above this variation set.
  Returntype : reference to list of Bio::EnsEMBL::Variation::VariationSet
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_super_VariationSets {
  my $self = shift;
  my $only_immediate = shift;
  
  # A database adaptor must be attached to this object   
  if(!$self->adaptor()) {
    warning('Cannot get super variation sets without attached adaptor');
    return [];
  }
  
  return $self->adaptor->fetch_all_by_sub_VariationSet($self,$only_immediate);
}

=head2 get_all_Variations

  Example    : print $vs->get_all_Variations();
  Description: Gets all variations belonging to this variation set and all of its subsets.
  Returntype : reference to list of Bio::EnsEMBL::Variation::Variation
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_Variations {
  my $self = shift;
  
  # A database adaptor must be attached to this object   
  if(!$self->adaptor()) {
    warning('Cannot get variations without attached adaptor');
    return [];
  }
  
  # Call the method in VariationAdaptor that will handle this
  my $variation_adaptor = $self->adaptor->db->get_VariationAdaptor();
  if(!$variation_adaptor) {
    warning('Could not get variation adaptor from database');
    return [];
  }

  # Get all variations from this set (and its subsets)
  return $variation_adaptor->fetch_all_by_VariationSet($self);
}


=head2 get_all_StructuralVariations

  Example    : print $vs->get_all_StructuralVariations();
  Description: Gets all structural variations belonging to this variation set and all of its subsets.
  Returntype : reference to list of Bio::EnsEMBL::Variation::StructuralVariation
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_StructuralVariations {
  my $self = shift;
  
  # A database adaptor must be attached to this object   
  if(!$self->adaptor()) {
    warning('Cannot get structural variations without attached adaptor');
    return [];
  }
  
  # Call the method in StructuralVariationAdaptor that will handle this
  my $sv_adaptor = $self->adaptor->db->get_StructuralVariationAdaptor();
  if(!$sv_adaptor) {
    warning('Could not get structural variation adaptor from database');
    return [];
  }

  # Get all variations from this set (and its subsets)
  return $sv_adaptor->fetch_all_by_VariationSet($self);
}


=head2 get_Variation_Iterator

  Example    : my $var_iterator = $vs->get_Variation_Iterator;
  Description: Gets an iterator over  all variations belonging to this 
	             variation set and all of its subsets.
  Returntype : Bio::EnsEMBL::Utils::Iterator
  Exceptions : none
  Caller     : general
  Status     : Experimental

=cut

sub get_Variation_Iterator {
    my $self = shift;
  
    # A database adaptor must be attached to this object   
    unless ($self->adaptor) {
        warning('Cannot get variations without attached adaptor');
        return Bio::EnsEMBL::Utils::Iterator->new;
    }
  
    # Call the method in VariationAdaptor that will handle this
    my $variation_adaptor = $self->adaptor->db->get_VariationAdaptor();
    
    unless ($variation_adaptor) {
        warning('Could not get variation adaptor from database');
        return Bio::EnsEMBL::Utils::Iterator->new;
    }

    # Get an iterator over variations from this set (and its subsets)
    return $variation_adaptor->fetch_Iterator_by_VariationSet($self);
}


=head2 get_StructuralVariation_Iterator

  Example    : my $sv_iterator = $vs->get_StructuralVariation_Iterator;
  Description: Gets an iterator over all structural variations belonging to this 
	             variation set and all of its subsets.
  Returntype : Bio::EnsEMBL::Utils::Iterator
  Exceptions : none
  Caller     : general
  Status     : Experimental

=cut

sub get_StructuralVariation_Iterator {
    my $self = shift;
  
    # A database adaptor must be attached to this object   
    unless ($self->adaptor) {
        warning('Cannot get structural variations without attached adaptor');
        return Bio::EnsEMBL::Utils::Iterator->new;
    }
  
    # Call the method in StructuralVariationAdaptor that will handle this
    my $sv_adaptor = $self->adaptor->db->get_StructuralVariationAdaptor();
    
    unless ($sv_adaptor) {
        warning('Could not get dructural variation adaptor from database');
        return Bio::EnsEMBL::Utils::Iterator->new;
    }

    # Get an iterator over Structural variations from this set (and its subsets)
    return $sv_adaptor->fetch_Iterator_by_VariationSet($self);
}


=head2 get_all_VariationFeatures_by_Slice

  Arg [1]    : Bio::EnsEMBL:Variation::Slice $slice
  Example    : my @vfs =
@{$vs->get_all_VariationFeatures_by_Slice($slice)};
  Description: Retrieves all variation features in a slice that belong to 
			   this set.
  Returntype : reference to list Bio::EnsEMBL::Variation::VariationFeature
  Exceptions : throw on bad argument
  Caller     : general
  Status     : Stable

=cut

sub get_all_VariationFeatures_by_Slice {
  my $self = shift;
  my $slice = shift;
  
  if(!$self->adaptor()) {
    warning('Cannot get variation features without attached adaptor');
    return [];
  }
  
  my $vfa = $self->adaptor->db->get_VariationFeatureAdaptor();
  if(!$vfa) {
    warning('Could not get variation feature adaptor from database');
    return [];
  }
  return $vfa->fetch_all_by_Slice_VariationSet($slice,$self);
}

=head2 name

  Arg [1]    : string $name
  Example    : print $vs->name();
  Description: Getter/Setter for the name of this VariationSet
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub name {
  my $self = shift;
  my $name = shift;
  
  $self->{'name'} = $name if (defined($name));
  
  return $self->{'name'};
}

=head2 short_name

  Arg [1]    : string $short_name
  Example    : print $vs->short_name();
  Description: Getter/Setter for the short name of this VariationSet
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub short_name {
  my $self = shift;
  my $short_name = shift;
  
  $self->{'short_name'} = $short_name if (defined($short_name));
  
  return $self->{'short_name'};
}

# API-internal subroutine to get the bitvalue of this set's id and all of its subsets (unless specifically indicated not to)
sub _get_bitvalue {
  my $self = shift;
  my @args = @_;
  
  # If the subsets should be exluded, call the subroutine in the adaptor and return the result. No caching.
  if (@args) {
    return $self->adaptor->_get_bitvalue($self,@args);
  }
  
  # Check if we have cached the bitvalue (including subsets), otherwise get it and cache it
  unless (exists($self->{'_bitvalue'})) {
    $self->{'_bitvalue'} = $self->adaptor->_get_bitvalue($self);
  }
  
  # Return the cached value
  return $self->{'_bitvalue'};
}

1;
