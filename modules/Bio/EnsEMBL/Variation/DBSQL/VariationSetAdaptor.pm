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

#
# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::VariationSetAdaptor
#
#

=head1 NAME

Bio::EnsEMBL::Variation::DBSQL::VariationSetAdaptor

=head1 SYNOPSIS

  $db = Bio::EnsEMBL::Variation::DBSQL::DBAdaptor->new(...);

  $vsa = $db->get_VariationSetAdaptor();

  # retrieve a variation set by its name
  $vs = $vsa->fetch_by_name('Phenotype-associated variations');

  # retrieve a variation set by its internal identifier
  $vs = $vsa->fetch_by_dbID(12);

  # retrieve all variation sets which a variation is a part of
  @vs = @{$vsa->fetch_all_by_Variation($var)};


=head1 DESCRIPTION

This adaptor provides database connectivity for VariationSet objects.
VariationSets may be retrieved from the Ensembl variation database by
several means using this module.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::VariationSetAdaptor;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Scalar qw(assert_ref wrap_array);

use Bio::EnsEMBL::Variation::VariationSet;

our @ISA = ('Bio::EnsEMBL::DBSQL::BaseAdaptor');

our $MAX_VARIATION_SET_ID = 64;

=head2 fetch_all_top_VariationSets

  Example    : $vs = $vs_adaptor->fetch_all_top_VariationSets();
  Description: Retrieves all VariationSet objects that are 'toplevel', 
               i.e. they are not subsets of any other variation set.
  Returntype : istref of Bio::EnsEMBL::Variation::VariationSet
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_top_VariationSets {
    my $self = shift;

    #Add a constraint to only get the sets that don't have any parent sets
    my $constraint = qq{
        NOT EXISTS (
            SELECT
                *
            FROM
                variation_set_structure vss
            WHERE
                vss.variation_set_sub = vs.variation_set_id
        )
    };
    
    #Get the results from generic fetch method
    return $self->generic_fetch($constraint);
    
}

=head2 fetch_all_by_sub_VariationSet

  Arg [1]    : Bio::EnsEMBL::Variation::VariationSet $sub
  Arg [2]    : (optional) boolean $only_immediate
               If true, only the direct supersets of this variation set will be fetched. The default behaviour is
               to recursively fetch all supersets.
  Example    : @vs_supersets = @{$vs_adaptor->fetch_all_by_sub_VariationSet($vs)};
  Description: Retrieves all VariationSets that are direct supersets of the specified VariationSet.
  Returntype : listref of Bio::EnsEMBL::Variation::VariationSet
  Exceptions : throw if sub arg is not valid
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_sub_VariationSet {
  my $self = shift;
  my $set = shift;
  my $only_immediate = shift;

  #Check the input set
  assert_ref($set,'Bio::EnsEMBL::Variation::VariationSet');
  
# First, get all VariationSets that are direct supersets of this one

  my $dbID = $set->dbID();
  
  my $stmt = qq{
    SELECT
      vss.variation_set_super
    FROM
      variation_set_structure vss
    WHERE
      vss.variation_set_sub = ?
  };
  my $sth = $self->prepare($stmt);
  $sth->execute($dbID);

  my %vs;
  while (my $result = $sth->fetchrow_arrayref()) {
# For each superset, fetch all of its supersets, unless specifically told not to
    my $vs_sup = $self->fetch_by_dbID($result->[0]);
    $vs{$vs_sup->dbID()} = $vs_sup;
    if (!defined($only_immediate)) {
      foreach my $v (@{$self->fetch_all_by_sub_VariationSet($vs_sup)}) {
        $vs{$v->dbID()} = $v;
      }
    }
  }
  
  my @res = values(%vs);
  
  return \@res;
}

=head2 fetch_all_by_super_VariationSet

  Arg [1]    : Bio::EnsEMBL::Variation::VariationSet $super
  Arg [2]    : (optional) boolean $only_immediate
               If true, only the direct subsets of this variation set will be fetched. The default behaviour is
               to recursively fetch all subsets.
  Example    : @vs_subsets = @{$vs_adaptor->fetch_all_by_super_VariationSet($vs)};
  Description: Retrieves all VariationSets that are subsets of the specified VariationSet.
  Returntype : listref of Bio::EnsEMBL::Variation::VariationSet
  Exceptions : throw if super arg is not valid
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_super_VariationSet {
  my $self = shift;
  my $set = shift;
  my $only_immediate = shift;
  
  #Check the input set
  assert_ref($set,'Bio::EnsEMBL::Variation::VariationSet');
  
# First, get all VariationSets that are direct subsets of this one

  my $dbID = $set->dbID();
  
  my $stmt = qq{
    SELECT
      vss.variation_set_sub
    FROM
      variation_set_structure vss
    WHERE
      vss.variation_set_super = ?
  };
  my $sth = $self->prepare($stmt);
  $sth->execute($dbID);

  my %vs;
  while (my $result = $sth->fetchrow_arrayref()) {
# For each subset, fetch all of its subsets unless specifically told not to
    my $vs_sub = $self->fetch_by_dbID($result->[0]);
    $vs{$vs_sub->dbID()} = $vs_sub;
    if (!defined($only_immediate)) {
      foreach my $v (@{$self->fetch_all_by_super_VariationSet($vs_sub)}) {
        $vs{$v->dbID()} = $v;
      }
    }
  }
  
  my @res = values(%vs);
  
  return \@res;
}


=head2 fetch_by_name

  Arg [1]    : string $name
  Example    : $vg = $vga->fetch_by_name('Phenotype-associated variations');
  Description: Retrieves a variation set by its name.
  Returntype : Bio::EnsEMBL::Variation::VariationSet
  Exceptions : throw if name argument is not provided
  Caller     : general
  Status     : Stable

=cut

sub fetch_by_name {
    my $self = shift;
    my $name = shift;

    throw('name argument expected') unless (defined($name));

    # Add a constraint on the name column and bind the name to it
    my $constraint = qq{ vs.name LIKE ? };
    $self->bind_param_generic_fetch($name,SQL_VARCHAR);
    
    #Call the generic fetch method
    my $result = wrap_array($self->generic_fetch($constraint));
    
    # Return the result
    return undef unless (scalar(@{$result}));
    return $result->[0];
}

=head2 fetch_by_short_name

  Arg [1]    : string $name
  Example    : $vg = $vga->fetch_by_short_name('ph_variants');
  Description: Retrieves a variation set by its short name.
  Returntype : Bio::EnsEMBL::Variation::VariationSet
  Exceptions : throw if short name argument is not provided
  Caller     : general

=cut

sub fetch_by_short_name {
    my $self = shift;
    my $name = shift;

    throw('short name argument expected') unless (defined($name));

    #Get the attrib_id corresponding to the 'short_name' type and specified name
    my $aa = $self->db->get_AttributeAdaptor();
    my $attrib_id = $aa->attrib_id_for_type_value($self->_short_name_attrib_type_code(),$name);
    return undef unless (defined($attrib_id));
    
    # Add a constraint on the short_name_attrib_id column and bind the name to it
    my $constraint = qq{ vs.short_name_attrib_id = ? };
    $self->bind_param_generic_fetch($attrib_id,SQL_INTEGER);
    
    #Call the generic fetch method
    my $result = wrap_array($self->generic_fetch($constraint));
    
    # Return the result
    return undef unless (scalar(@{$result}));
    return $result->[0];
}


=head2 fetch_all_by_Variation

  Arg [1]    : Bio::EnsEMBL::Variation::Variation
  Example    : my $vgs = $vga->fetch_all_by_Variation($var);
  Description: Retrieves all variation sets which a particular variation
               is present in.
  Returntype : reference to list of Bio::EnsEMBL::Variation::VariationSets
  Exceptions : throw on incorrect argument
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_Variation {
    my $self = shift;
    my $var  = shift;

    assert_ref($var,'Bio::EnsEMBL::Variation::Variation');

  my $cols = join(',',$self->_columns());
  my $stmt = qq{
    SELECT
      $cols
    FROM
      variation_set vs,
      variation_set_variation vsv
    WHERE
      vs.variation_set_id = vsv.variation_set_id AND
      vsv.variation_id = ?
  };

  my $sth = $self->prepare($stmt);
  $sth->bind_param(1,$var->dbID,SQL_INTEGER);
  $sth->execute();

  my $result = $self->_objs_from_sth($sth);
  $sth->finish();

# Fetch all supersets of the returned sets as well. Since a variation may occur at several places in a hierarchy
# which will cause duplicated data, store variation sets in a hash with dbID as key.
  my %sets;
  foreach my $set (@{$result}) {
    $sets{$set->dbID()} = $set;
    foreach my $sup (@{$self->fetch_all_by_sub_VariationSet($set)}) {
      $sets{$sup->dbID()} = $sup;
    }
  }
  
  my @res = values %sets;
  return \@res;
}



=head2 fetch_all_by_Variation_super_VariationSet

  Arg [1]    : Bio::EnsEMBL::Variation::Variation
  Arg [2]    : Bio::EnsEMBL::Variation::VariationSet
  Example    : my $vs = $vsa->fetch_all_by_Variation_super_VariationSet($var, $var_set);
  Description: Retrieves all variation sets which a particular variation
               is present in, which are subsets of a given super set.
  Returntype : reference to list of Bio::EnsEMBL::Variation::VariationSets
  Exceptions : throw on incorrect argument
  Caller     : general
  Status     : Experimental

=cut

sub fetch_all_by_Variation_super_VariationSet {

    my $self      = shift;
    my $var       = shift;
    my $superset  = shift;

    assert_ref($var,'Bio::EnsEMBL::Variation::Variation');
    assert_ref($superset,'Bio::EnsEMBL::Variation::VariationSet');

  my $cols = join(',',$self->_columns());
  my $stmt = qq{
    SELECT
      $cols
    FROM
      variation_set vs,
      variation_set vssup,
      variation_set_structure vss,
      variation_set_variation vsv
    WHERE
      vss.variation_set_super = ? AND 
      vss.variation_set_sub =  vs.variation_set_id AND
      vs.variation_set_id = vsv.variation_set_id AND
      vsv.variation_id = ?
  };

  my $sth = $self->prepare($stmt);

  $sth->execute($superset->dbID(), $var->dbID() );

  my $result = $self->_objs_from_sth($sth);
  $sth->finish();

# store variation sets in a hash with dbID as key to remove duplicates.
  my %sets;
  foreach my $set (@{$result}) {
    $sets{$set->dbID()} = $set;   
  }
  
  my @res = values %sets;
  return \@res;
}


=head2 fetch_all_by_StructuralVariation

  Arg [1]    : Bio::EnsEMBL::Variation::StructuralVariation
  Example    : my $vss = $vsa->fetch_all_by_StructuralVariation($sv);
  Description: Retrieves all variation sets which a particular structural variation
               is present in.
  Returntype : reference to list of Bio::EnsEMBL::Variation::VariationSets
  Exceptions : throw on incorrect argument
  Caller     : general
  Status     : Stable

=cut

sub fetch_all_by_StructuralVariation {
    my $self = shift;
    my $var  = shift;

    assert_ref($var,'Bio::EnsEMBL::Variation::StructuralVariation');

  my $cols = join(',',$self->_columns());
  my $stmt = qq{
    SELECT
      $cols
    FROM
      variation_set vs,
      variation_set_structural_variation vssv
    WHERE
      vs.variation_set_id = vssv.variation_set_id AND
      vssv.structural_variation_id = ?
  };

  my $sth = $self->prepare($stmt);
  $sth->bind_param(1,$var->dbID,SQL_INTEGER);
  $sth->execute();

  my $result = $self->_objs_from_sth($sth);
  $sth->finish();

# Fetch all supersets of the returned sets as well. Since a variation may occur at several places in a hierarchy
# which will cause duplicated data, store variation sets in a hash with dbID as key.
  my %sets;
  foreach my $set (@{$result}) {
    $sets{$set->dbID()} = $set;
    foreach my $sup (@{$self->fetch_all_by_sub_VariationSet($set)}) {
      $sets{$sup->dbID()} = $sup;
    }
  }
  
  my @res = values %sets;
  return \@res;
}


# An API-internal subroutine for getting the bitvalue of the specified variation_set and (unless specifically indicated) its subsets
sub _get_bitvalue {
  my $self = shift;
  my $set = shift;
  return sprintf("power(2, %i)", $set->dbID() - 1);
}

# API-internal method for getting the attrib_type code used for short names
sub _short_name_attrib_type_code {
    return q{short_name};
}

sub _columns {
  return qw( vs.variation_set_id vs.name vs.description vs.short_name_attrib_id );
}
sub _tables {
  return ( ['variation_set','vs'] );
}
sub _default_where_clause {
  return '1';
}

sub _objs_from_sth {
  my $self = shift;
  my $sth  = shift;

  my ($vs_id, $name, $description, $short_name_attrib_id);
  $sth->bind_columns(\$vs_id, \$name, \$description, \$short_name_attrib_id);

  my @results;
  my ($cur_vs, $cur_vs_id);
    my $aa = $self->db->get_AttributeAdaptor();
    
# Construct all variation sets

  while($sth->fetch()) {
    if (!defined($cur_vs) || $vs_id != $cur_vs_id) {
      $cur_vs = Bio::EnsEMBL::Variation::VariationSet->new
        (
          -dbID => $vs_id,
          -adaptor => $self,
          -name    => $name,
          -description  => $description,
          -short_name   => $aa->attrib_value_for_id($short_name_attrib_id)
        );
      $cur_vs_id = $vs_id;
      push(@results,$cur_vs);
    }
  }

  return \@results;
}


1;
