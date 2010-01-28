# Ensembl module for Bio::EnsEMBL::Variation::VariationSet
#
# Copyright (c) 2010 Ensembl
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

=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::Variation::VariationSet;

use Bio::EnsEMBL::Storable;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw deprecate warning);

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

  Example    :
    $ag = Bio::EnsEMBL::Variation::VariationSet->new
      (
       -dbID => 12,
       -adaptor => $var_set_adaptor,
       -name   => 'Phenotype-associated variations',
       -description => 'Variations that have been associated with a phenotype'
      );
  Description: Constructor.  Instantiates a new VariationSet
  Returntype : Bio::EnsEMBL::Variation::VariationSet
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub new {
  my $class = shift;

  my ($dbID, $adaptor, $name, $description) =
    rearrange([qw(DBID ADAPTOR NAME DESCRIPTION)], @_);
  
  return bless {'dbID' => $dbID,
                'adaptor' => $adaptor,
                'name' => $name,
                'description' => $description}, $class;
}

=head2 description

  Arg [1]    : string $description
  Example    : print $vs->description();
  Description: Getter/Setter for the description of this VariationSet
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

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
  Status     : At Risk

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
  Status     : At Risk

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
  Status     : At Risk

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
  Status     : At Risk

=cut

sub name {
  my $self = shift;
  my $name = shift;
  
  $self->{'name'} = $name if (defined($name));
  
  return $self->{'name'};
}

1;
