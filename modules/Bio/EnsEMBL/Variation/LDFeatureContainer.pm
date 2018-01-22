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

# Ensembl module for Bio::EnsEMBL::Variation::LDFeatureContainer
#
#


=head1 NAME

Bio::EnsEMBL::Variation::LDFeatureContainer - A container with all the ld values for quick access

=head1 SYNOPSIS

  # LDFeature Container representing all the LD values for a certain contig
  $ldContainer = Bio::EnsEMBL::Variation::LDFeatureContainer->new(
    -name   => NT_085213,
    -ldContainer => $ldhash,
    -variation_features => $vfhash);

  # get the d_prime values for a certain pair of variation_features
  d_prime = $ldContainer->get_d_prime($vf1, $vf2);

  # get list variants in the container
  $variations = $ldContainer->get_variations();

=head1 DESCRIPTION

This is a class representing the LD information for a certain region
from the ensembl-variation database.
See B<Bio::EnsEMBL::Variation::Variation>.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::LDFeatureContainer;

use Bio::EnsEMBL::Storable;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument  qw(rearrange);

use vars qw(@ISA);
use Scalar::Util qw(looks_like_number);

@ISA = qw(Bio::EnsEMBL::Storable);

=head2 new

  Arg [-NAME] :
    string - name of the feature object that the LD container refers to. (chr1,NT_08542,...)

  Arg [-LDCONTAINER] :
    reference - hash containing all the LD information present, with the key
    (vf1-vf2) to access the information

  Arg [-VARIATIONFEATURES] :
    reference - hash containing all the Bio::EnsEMBL::Variation::VariationFeature objects that are present in the Container
    
  Example    :
  $ldContainer = Bio::EnsEMBL::Variation::LDFeatureContainer->new(
    -name => 'chr1'
    -ldContainer => {
      'vf1-vf2' => {
        'population_id_1' => {
          'd_prime' => 0.5,
          'r2'      => 0.421,
          'sample_count' => 120
        },
        'population_id_2' => {
          'd_prime' => 0.3,
          'r2'     => 0.321,
          'sample_count' => 35
        }
      }
    }
    -variationFeatures => hash of Bio::EnsEMBL::Variation::VariationFeature
  );


  Description: Constructor. Instantiates a new LDFeatureContainer object.
  Returntype : Bio::EnsEMBL::Variation::LDFeatureContainer
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my ($ldContainer,$name,$pos2vf,$pos2name,$slices,$adaptor) =  rearrange([qw(LDCONTAINER NAME POS2VF POS2NAME SLICES ADAPTOR)], @_);
  if (defined($ldContainer) && ref($ldContainer ne 'HASH')){
    throw("Reference to a hash object expected as a LDContainer");
  }
  $ldContainer ||= {};
  $pos2name ||= {};

  my $self = bless {
    'name' => $name,
    'ldContainer' => $ldContainer,
    'pos2name' => $pos2name,
    'slices' => $slices,
    'adaptor' => $adaptor,
  }, $class;

  # only add these keys if it they are properly populated
  # makes the lazy-load later easier
  $self->{pos2vf} = $pos2vf if $pos2vf && scalar keys %$pos2vf;

  return $self;
}


=head2 name

  Arg [1]    : string $newval (optional)
               The new value to set the name attribute to
  Example    : $name = $obj->name()
  Description: Getter/Setter for the name attribute.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub name{
  my $self = shift;
  return $self->{'name'} = shift if(@_);
  return $self->{'name'};
}


### THIS SHOULD REALLY BE NAMED get_all_Variations
### IS IT STILL USED???
=head2 get_variations

    Example     : $variations = $obj->get_variations()
    Description : Gets all the variation objects contained in the LDFeatureContainer
    ReturnType  : list of Bio::EnsEMBL::Variation::Variation
    Exceptions  : none
    Caller      : general
    Status      : At Risk

=cut

sub get_variations {
  my $self = shift;
  my $pos2vf = $self->_pos2vf();
  return [map {$pos2vf->{$_}->variation()} sort {$a <=> $b} keys %$pos2vf];
}


=head2 get_r_square

    Arg [1]     : Bio::EnsEMBL::Variation::VariationFeature $variationFeature
    Arg [2]     : Bio::EnsEMBL::Variation::VariationFeature $variationFeature
    Arg [3]     : (optional) int - population_id of population you want to get the r_square value
    Example     : $r_square = $obj->get_r_square($vf1,$vf2,$population_id);
    Description : Get the r_square value for a pair of variation features in the given population. If no population is provided,
    return the r_square for the default population with more sample counts (in case more than 1)
    ReturnType  : float
    Exceptions  : throw on incorrect arguments
    Caller      : general
    Status      : At Risk

=cut

sub get_r_square {
  my $self = shift;
  my $vf1 = shift;
  my $vf2 = shift;
  my $population_id = shift;

  $population_id ||= 0; #in case no population provided, to avoid warning in the hash
  my $r_square;
  my $key;

  my $vf1_key = $self->_get_vf_key($vf1);
  my $vf2_key = $self->_get_vf_key($vf2);

  #check if the default poppulation has been calculated, otherwise, find it
  if (! defined $self->{'_default_population'}){
    $self->{'_default_population'} = $self->_get_major_population;
  }
  #first of all, check that both arguments have been properly provided
  if (!ref($vf1) || !$vf1->isa('Bio::EnsEMBL::Variation::VariationFeature') || !ref($vf2) || !$vf2->isa('Bio::EnsEMBL::Variation::VariationFeature')){
    throw("Bio::EnsEMBL::Variation::VariationFeature arguments expected");
  }
  else {
    #check if the ldContainer does not contain pairwise information for the variation features provided
    if (!exists $self->{'ldContainer'}->{$vf1_key.'-'.$vf2_key} && !exists $self->{'ldContainer'}->{$vf2_key.'-'.$vf1_key}){
      warning("variation features have no pairwise ld information");
    } 
    else {
      #find out the key in the ld Hash: vf1 - vf2 or vf2 - vf1
      if (exists $self->{'ldContainer'}->{$vf1_key.'-'.$vf2_key}){
        $key = $vf1_key.'-'.$vf2_key;
      }
      else {
        $key = $vf2_key.'-'.$vf1_key;
      }
      #and finally, if population provided or the only population
      if (exists $self->{'ldContainer'}->{$key}->{$population_id}){
        $r_square = $self->{'ldContainer'}->{$key}->{$population_id}->{'r2'}
      }
      else {
        if (exists $self->{'ldContainer'}->{$key}->{$self->{'_default_population'}}){
          #there was no population provided, return the r_square for the default population
          $r_square = $self->{'ldContainer'}->{$key}->{$self->{'_default_population'}}->{'r2'};		    
        }
        else {
          warning("variation features have no pairwise ld information for default population ", $self->{'_default_population'});
        }
      }
    }
  }
  return $r_square;
}


=head2 get_d_prime

    Arg [1]     : Bio::EnsEMBL::Variation::VariationFeature $variationFeature
    Arg [2]     : Bio::EnsEMBL::Variation::VariationFeature $variationFeature
    Arg [3]     : (optional) int - population_id of population you want to get the d_prime value
    Example     : $d_prime = $obj->get_d_prime($vf1, $vf2, $population_id);
    Description : Get the d_prime value for a pair of variation features for a known or unknown population. In case of an unknown population, the default
poulation is used    
    ReturnType  : float
    Exceptions  : throw on incorrect arguments
    Caller      : general
    Status      : At Risk

=cut

sub get_d_prime {
  my $self = shift;
  my $vf1 = shift;
  my $vf2 = shift;
  my $population_id = shift;

  $population_id ||= 0; #in case no population provided, to avoid warning in the hash
  my $d_prime;
  my $key;

  my $vf1_key = $self->_get_vf_key($vf1);
  my $vf2_key = $self->_get_vf_key($vf2);

  if (! defined $self->{'_default_population'}){
    $self->{'_default_population'} = $self->_get_major_population;
  }
  #first of all, check that both arguments have been properly provided
  if (!ref($vf1) || !$vf1->isa('Bio::EnsEMBL::Variation::VariationFeature') || !ref($vf2) || !$vf2->isa('Bio::EnsEMBL::Variation::VariationFeature')){
    throw("Bio::EnsEMBL::Variation::VariationFeature arguments expected");
  }
  else {
    #check if the ldContainer does not contain pairwise information for the variation features provided
    if (!exists $self->{'ldContainer'}->{$vf1_key.'-'.$vf2_key} && !exists $self->{'ldContainer'}->{$vf2_key.'-'.$vf1_key}){
      warning("variation features have no pairwise ld information");
    } 
    else {
      #find out the key in the ld Hash: vf1 - vf2 or vf2 - vf1
      if (exists $self->{'ldContainer'}->{$vf1_key.'-'.$vf2_key}){
        $key = $vf1_key.'-'.$vf2_key;
      }
      else {
        $key = $vf2_key.'-'.$vf1_key;
      }
      #and finally, if population provided or the only population
      if (exists $self->{'ldContainer'}->{$key}->{$population_id}){
        $d_prime = $self->{'ldContainer'}->{$key}->{$population_id}->{'d_prime'};
      }
      else {
        if (exists $self->{'ldContainer'}->{$key}->{$self->{'_default_population'}}){
          #there was no population provided, return the r_square for the default population
          $d_prime = $self->{'ldContainer'}->{$key}->{$self->{'_default_population'}}->{'d_prime'};
        }
        else {
          warning("variation features have no pairwise ld information for default population ", $self->{'_default_population'});
        }
      }
    }    
  }
  return $d_prime;
}


=head2 get_all_ld_values

    Arg [1]     : bool $names_only - set to a true value to populate hashes with names only, not VariationFeature objects also
                  Defaults to fetching VariationFeature objects too
    Example     : $ld_values = $obj->get_all_ld_values();
    Description : Get all the information contained in the LDFeatureContainer object
    ReturnType  : reference to list of hashes [{variation1 => Bio::EnsEMBL::Variation::VariationFeature, variation2=>Bio::EnsEMBL::Variation::VariationFeature, d_prime=>d_prime, r2=>r2, sample_count=>sample_count, population_id=>population_id}]
    Exceptions  : no exceptions
    Caller      : general
    Status      : At Risk

=cut

sub get_all_ld_values {
  my $self = shift;
  my $names_only = shift;

  my @ld_values; #contains ALL the ld values in the container

  if (! defined $self->{'_default_population'}){
    $self->{'_default_population'} = $self->_get_major_population;
  }

  # get these hashes for looking up names and VF objects
  my $pos2name = $self->_pos2name();
  my $pos2vf = $self->_pos2vf() unless $names_only;
  my $vf_name = $self->{'_vf_name'};
  foreach my $key_ld (keys %{$self->{'ldContainer'}}) {
    # contains a single ld value in the container {variation_feature variation_feature d_prime r2}
    my %ld_value;  

    # get the variation_features positions
    my ($vf1_pos, $vf2_pos) =  split /-/,$key_ld;

    if (exists $self->{'ldContainer'}->{$key_ld}->{$self->{'_default_population'}}) {
      
      # add the information to the ld_value hash
      if (!$names_only) {
        if ($vf_name) {
          if ($pos2vf->{$vf1_pos} && $pos2vf->{$vf1_pos}->{variation_name} eq $vf_name) {
            $ld_value{'variation1'} = $pos2vf->{$vf1_pos};
            $ld_value{'variation2'} = $pos2vf->{$vf2_pos};
          } else {
            $ld_value{'variation2'} = $pos2vf->{$vf1_pos};
            $ld_value{'variation1'} = $pos2vf->{$vf2_pos};
          }
        } else {
          $ld_value{'variation1'} = $pos2vf->{$vf1_pos};
          $ld_value{'variation2'} = $pos2vf->{$vf2_pos};
        }
        # $DB::single = 1 unless $ld_value{'variation1'} && $ld_value{'variation2'};
        next unless $ld_value{'variation1'} && $ld_value{'variation2'};
      }
      if ($vf_name) {
        if ($pos2name->{$vf1_pos} eq $vf_name) {
          $ld_value{'variation_name1'} = $pos2name->{$vf1_pos};
          $ld_value{'variation_name2'} = $pos2name->{$vf2_pos};
        } else {
          $ld_value{'variation_name2'} = $pos2name->{$vf1_pos};
          $ld_value{'variation_name1'} = $pos2name->{$vf2_pos};
        }  
      } else {
        $ld_value{'variation_name1'} = $pos2name->{$vf1_pos};
        $ld_value{'variation_name2'} = $pos2name->{$vf2_pos};
      }

      next unless $ld_value{'variation_name1'} && $ld_value{'variation_name2'};

      $ld_value{'d_prime'} = $self->{'ldContainer'}->{$key_ld}->{$self->{'_default_population'}}->{'d_prime'};
      $ld_value{'r2'} = $self->{'ldContainer'}->{$key_ld}->{$self->{'_default_population'}}->{'r2'};
      $ld_value{'sample_count'} = $self->{'ldContainer'}->{$key_ld}->{$self->{'_default_population'}}->{'sample_count'};
      $ld_value{'population_id'} = $self->{'_default_population'};
      push @ld_values, \%ld_value;
    }
  }

  return \@ld_values;
}


=head2 get_all_r_square_values

    Arg [1]     : bool $names_only - set to a true value to populate hashes with names only, not VariationFeature objects also
                  Defaults to fetching VariationFeature objects too
    Example     : $r_square_values = $obj->get_all_r_square_values();
    Description : Get all r_square values contained in the LDFeatureContainer object
    ReturnType  : reference to list of [{variation1=>Bio::EnsEMBL::Variation::VariationFeature, variation2=>Bio::EnsEMBL::Variation::VariationFeature, r2=>r2, population_id=>population_id}]
    Exceptions  : no exceptions
    Caller      : general
    Status      : At Risk

=cut

sub get_all_r_square_values {
  my $self = shift;
  return [map {delete $_->{d_prime}; $_} @{$self->get_all_ld_values(@_)}];
}


=head2 get_all_d_prime_values

    Arg [1]     : bool $names_only - set to a true value to populate hashes with names only, not VariationFeature objects also
                  Defaults to fetching VariationFeature objects too
    Example     : $d_prime_values = $obj->get_all_d_prime_values();
    Description : Get all d_prime values contained in the LDFeatureContainer object
    ReturnType  : reference to list of [{variation1=>Bio::EnsEMBL::Variation::VariationFeature, variation2=>Bio::EnsEMBL::Variation::VariationFeature, d_prime=>d_prime, population_id=>population_id}]
    Exceptions  : no exceptions
    Caller      : general
    Status      : At Risk

=cut

sub get_all_d_prime_values {
  my $self = shift;
  return [map {delete $_->{r2}; $_} @{$self->get_all_ld_values(@_)}];
}


=head2 get_all_populations

    Arg [1]     : (optional) Bio::EnsEMBL::Variation::VariationFeature $variationFeature
    Arg [2]     : (optional) Bio::EnsEMBL::Variation::VariationFeature $variationFeature
    Example     : $populations = $obj->get_all_populations($vf1,$vf2);
    Description : If no arguments provided, returns ALL the populations present in the container. When 2 variation features provided, returns the 
    population/populations where these variation features occurs
    ReturnType  : reference to list of int
    Exceptions  : throw on incorrect arguments
    Caller      : general
    Status      : At Risk

=cut

sub get_all_populations {
  my $self = shift;
  my $vf1 = shift;
  my $vf2 = shift;
  my %populations;
  my @populations;    
  my $key;

  #if no variation provided, return ALL the populations in the container
  if (! defined($vf1) && ! defined($vf2)) {
    foreach my $key (keys %{$self->{'ldContainer'}}){
      map {$populations{$_}++} keys %{$self->{'ldContainer'}->{$key}};
    }
    @populations = keys %populations;
  }
  else {
    #first, check if both arguments have been properly provided
    if (!ref($vf1) || !$vf1->isa('Bio::EnsEMBL::Variation::VariationFeature') || !ref($vf2) || !$vf2->isa('Bio::EnsEMBL::Variation::VariationFeature')){
      throw("Bio::EnsEMBL::Variation::VariationFeature arguments expected");
    }
    #if the variation_features are correct, return the list of populations
    else {
      my $vf1_key = $self->_get_vf_key($vf1);
      my $vf2_key = $self->_get_vf_key($vf2);

      #find out the key in the ld Hash: vf1 - vf2 or vf2 - vf1
      if (exists $self->{'ldContainer'}->{$vf1_key.'-'.$vf2_key}){
        $key = $vf1_key.'-'.$vf2_key;
      }
      else {
        $key = $vf2_key.'-'.$vf1_key;
      }	    
      @populations = keys %{$self->{'ldContainer'}->{$key}};
    }
  }

  return \@populations;
}

#returns from the container the population_id with the maximum number of pairwise_ld 
sub _get_populations {
  my $self = shift;
  my %populations;

  foreach my $key (keys %{$self->{'ldContainer'}}){
    map {$populations{$_}++} keys %{$self->{'ldContainer'}->{$key}};
  }
  my @sorted_populations = sort{$populations{$b} <=> $populations{$a}} keys %populations;
  return @sorted_populations;
}

sub _get_major_population { 
  my( $pop ) = $_[0]->_get_populations;
  return $pop;
}

## gets the slices associated with this container
sub _slices {
  return $_[0]->{slices};
}

## gets the hashref mapping position to variant name
sub _pos2name {
  return $_[0]->{pos2name};
}

## gets the hashref mapping position to VariationFeature object
## lazy-loaded using slices if not populated on object creation
sub _pos2vf {
  my $self = shift;
  if (!$self->adaptor) {
    throw('Adaptor is not available. Cannot fetch VariationFeatures.');
  }
  my $vf_adaptor = $self->adaptor->db->get_VariationFeatureAdaptor();
  if(!exists($self->{pos2vf})) {
    
    # dont bother if there's nothing in the container
    return $self->{pos2vf} = {} unless keys %{$self->{ldContainer}};

    my %pos2vf = ();
    my @slice_objects = ();
    foreach my $slice (@{$self->_slices}) {
      if (ref($slice) eq 'ARRAY') {
        foreach (@$slice) {
          push @slice_objects, $_;
        }
      } else {
        push @slice_objects, $slice;
      }
    }
    foreach my $slice (@slice_objects) {
      my $vfs = $vf_adaptor->fetch_all_by_Slice($slice);
      $_->display_consequence for grep {!looks_like_number($_->dbID)} @$vfs;
      my $region_Slice = $slice->seq_region_Slice();

      my $pos2name = $self->_pos2name();

      my %names = map {$_ => 1} values %$pos2name;
      my @filtered_vfs = ();
      foreach my $vf (@$vfs) {
        if ($names{$vf->variation_name}) {
          push @filtered_vfs, $vf;
        }
      }
      # the sort here "favours" variants from dbSNP since it is assumed dbSNP will have source_id = 1
      # otherwise we'd get co-located COSMIC IDs overwriting the dbSNP ones

      $pos2vf{$_->[1]} = $_->[0] for
        grep {$pos2name->{$_->[1]}}
        map {[$_, $_->start]}
        map {$_->transfer($region_Slice)}
        sort {($b->{_source_id} || 1) <=> ($a->{_source_id} || 1)}
        @filtered_vfs;

    }

    $self->{pos2vf} = \%pos2vf;
  }

  return $self->{pos2vf};
}

## this is a helper method used when the VF does not have a slice attached
## mainly to allow the test VFs w/o slices to work OK
sub _get_vf_key {
  my $self = shift;
  my $vf = shift;
  return $vf->{slice} ? $vf->seq_region_start : $vf->start;
}

1;

