=head1 LICENSE

 Copyright (c) 1999-2013 The European Bioinformatics Institute and
 Genome Research Limited.  All rights reserved.

 This software is distributed under a modified Apache license.
 For license details, please see

   http://www.ensembl.org/info/about/legal/code_licence.html

=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <dev@ensembl.org>.

 Questions may also be sent to the Ensembl help desk at
 <helpdesk@ensembl.org>.

=cut

# Ensembl module for Bio::EnsEMBL::Variation::LDFeatureContainer
#
# Copyright (c) 2004 Ensembl
#


=head1 NAME

Bio::EnsEMBL::Variation::LDFeatureContainer - A container with all the ld values for quick access

=head1 SYNOPSIS

    # LDFeature Container representing all the LD values for a certain contig
    $ldContainer = Bio::EnsEMBL::Variation::LDFeatureContainer->new
       (-name   => NT_085213,
        -ldContainer => $ldhash,
        -variation_features => $vfhash);

    ...

    #get the d_prime values for a certain pair of variation_features
    d_prime = $ldContainer->get_d_prime($variation_feature_1,$variation_feature_2);
    #get the list of variation in the container
    $variations = $ldContainer->get_variations();

    ...

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

@ISA = qw(Bio::EnsEMBL::Storable);

=head2 new

  Arg [-NAME] :
    string - name of the feature object that the LD container refers to. (chr1,NT_08542,...)

  Arg [-LDCONTAINER] :
    reference - hash containing all the LD information present, with the key
    (variation_feature_1-variation_feature_2) to access the information

  Arg [-VARIATIONFEATURES] :
    reference - hash containing all the Bio::EnsEMBL::Variation::VariationFeature objects that are present in the Container
    
  Example    :
    $ldContainer = Bio::EnsEMBL::Variation::LDFeatureContainer->new(
	  -name => 'chr1'
      -ldContainer => {
		'variation_feature_1-variation_feature_2' => {
		  'sample_id_1' => {
			'd_prime' => 0.5,
			'r2'      => 0.421,
			'sample_count' => 120
		  },
		  'sample_id_2' => {
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

  my ($ldContainer,$name,$variationFeatures) =  rearrange([qw(LDCONTAINER NAME VARIATIONFEATURES)], @_);
  if (defined($ldContainer) && ref($ldContainer ne 'HASH')){
      throw("Reference to a hash object expected as a LDContainer");
  }
  $ldContainer ||= {};

  return bless {'name' => $name,
		'ldContainer' => $ldContainer,
		'variationFeatures' => $variationFeatures}, $class;
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


=head2 get_variations

    Example     : $variations = $obj->get_variations()
    Description : Gets all the variation objects contained in the LDFeatureContainer
    ReturnType  : list of Bio::EnsEMBL::Variation::Variation
    Exceptions  : none
    Caller      : general
    Status      : At Risk

=cut

sub get_variations{
    my $self = shift;
    my @variations;

    foreach my $variation_feature (keys %{$self->{'variationFeatures'}}){
	push @variations,$self->{'variationFeatures'}->{$variation_feature}->variation();
    }    
    return \@variations;
}

=head2 get_r_square

    Arg [1]     : Bio::EnsEMBL::Variation::VariationFeature $variationFeature
    Arg [2]     : Bio::EnsEMBL::Variation::VariationFeature $variationFeature
    Arg [3]     : (optional) int - sample_id of population you want to get the r_square value
    Example     : $r_square = $obj->get_r_square($vf1,$vf2,$sample_id);
    Description : Get the r_square value for a pair of variation features in the given population. If no population is provided,
    return the r_square for the default population with more sample counts (in case more than 1)
    ReturnType  : float
    Exceptions  : throw on incorrect arguments
    Caller      : general
    Status      : At Risk

=cut

sub get_r_square{
    my $self = shift;
    my $variation_feature_1 = shift;
    my $variation_feature_2 = shift;
    my $sample_id = shift;
    
    $sample_id ||= 0; #in case no population provided, to avoid warning in the hash
    my $r_square;
    my $key;

    #check if the default poppulation has been calculated, otherwise, find it
    if (! defined $self->{'_default_population'}){
	$self->{'_default_population'} = $self->_get_major_population;
    }
    #first of all, check that both arguments have been properly provided
    if (!ref($variation_feature_1) || !$variation_feature_1->isa('Bio::EnsEMBL::Variation::VariationFeature') || !ref($variation_feature_2) || !$variation_feature_2->isa('Bio::EnsEMBL::Variation::VariationFeature')){
	throw("Bio::EnsEMBL::Variation::VariationFeature arguments expected");
    }
    else{
	#check if the ldContainer does not contain pairwise information for the variation features provided
	if (!exists $self->{'ldContainer'}->{$variation_feature_1->dbID() . '-' . $variation_feature_2->dbID()} && !exists $self->{'ldContainer'}->{$variation_feature_2->dbID() . '-' . $variation_feature_1->dbID()}){
	    warning("variation features have no pairwise ld information");
	} 
	else{
	    #find out the key in the ld Hash: vf1 - vf2 or vf2 - vf1
	    if (exists $self->{'ldContainer'}->{$variation_feature_1->dbID() . '-' . $variation_feature_2->dbID()}){
		$key = $variation_feature_1->dbID() . '-' . $variation_feature_2->dbID();
	    }
	    else{
		$key = $variation_feature_2->dbID() . '-' . $variation_feature_1->dbID();
	    }
	    #and finally, if population provided or the only population
	    if (exists $self->{'ldContainer'}->{$key}->{$sample_id}){
		$r_square = $self->{'ldContainer'}->{$key}->{$sample_id}->{'r2'}
	    }
	    else{
		if (exists $self->{'ldContainer'}->{$key}->{$self->{'_default_population'}}){
		    #there was no population provided, return the r_square for the default population
		    $r_square = $self->{'ldContainer'}->{$key}->{$self->{'_default_population'}}->{'r2'};		    
		}
		else{
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
    Arg [3]     : (optional) int - sample_id of population you want to get the d_prime value
    Example     : $d_prime = $obj->get_d_prime($vf1,$vf2,$sample_id);
    Description : Get the d_prime value for a pair of variation features for a known or unknown population. In case of an unknown population, the default
poulation is used    
    ReturnType  : float
    Exceptions  : throw on incorrect arguments
    Caller      : general
    Status      : At Risk

=cut

sub get_d_prime{
    my $self = shift;
    my $variation_feature_1 = shift;
    my $variation_feature_2 = shift;
    my $sample_id = shift;

    $sample_id ||= 0; #in case no population provided, to avoid warning in the hash
    my $d_prime;
    my $key;

    if (! defined $self->{'_default_population'}){
	$self->{'_default_population'} = $self->_get_major_population;
    }
    #first of all, check that both arguments have been properly provided
    if (!ref($variation_feature_1) || !$variation_feature_1->isa('Bio::EnsEMBL::Variation::VariationFeature') || !ref($variation_feature_2) || !$variation_feature_2->isa('Bio::EnsEMBL::Variation::VariationFeature')){
	throw("Bio::EnsEMBL::Variation::VariationFeature arguments expected");
    }
    else{
	#check if the ldContainer does not contain pairwise information for the variation features provided
	if (!exists $self->{'ldContainer'}->{$variation_feature_1->dbID() . '-' . $variation_feature_2->dbID()} && !exists $self->{'ldContainer'}->{$variation_feature_2->dbID() . '-' . $variation_feature_1->dbID()}){
	    warning("variation features have no pairwise ld information");
	} 
	else{
	    #find out the key in the ld Hash: vf1 - vf2 or vf2 - vf1
	    if (exists $self->{'ldContainer'}->{$variation_feature_1->dbID() . '-' . $variation_feature_2->dbID()}){
		$key = $variation_feature_1->dbID() . '-' . $variation_feature_2->dbID();
	    }
	    else{
		$key = $variation_feature_2->dbID() . '-' . $variation_feature_1->dbID();
	    }
	    #and finally, if population provided or the only population
	    if (exists $self->{'ldContainer'}->{$key}->{$sample_id}){
		$d_prime = $self->{'ldContainer'}->{$key}->{$sample_id}->{'d_prime'};
	    }
	    else{
		if (exists $self->{'ldContainer'}->{$key}->{$self->{'_default_population'}}){
		    #there was no population provided, return the r_square for the default population
		    $d_prime = $self->{'ldContainer'}->{$key}->{$self->{'_default_population'}}->{'d_prime'};
		}
		else{
		    warning("variation features have no pairwise ld information for default population ", $self->{'_default_population'});
		}
	    }
	}	    
    }
    return $d_prime;
}


=head2 get_all_ld_values

    Example     : $ld_values = $obj->get_all_ld_values();
    Description : Get all the information contained in the LDFeatureContainer object
    ReturnType  : reference to list of hashes [{variation1 => Bio::EnsEMBL::Variation::VariationFeature, variation2=>Bio::EnsEMBL::Variation::VariationFeature, d_prime=>d_prime, r2=>r2, sample_count=>sample_count, sample_id=>population_sample_id}]
    Exceptions  : no exceptions
    Caller      : general
    Status      : At Risk

=cut


sub get_all_ld_values{
    my $self = shift;
    my @ld_values; #contains ALL the ld values in the container

    #the keys in the ldContainer hash
    my $variation_feature_id_1;
    my $variation_feature_id_2;
    
    if (! defined $self->{'_default_population'}){
	$self->{'_default_population'} = $self->_get_major_population;
    }
    foreach my $key_ld (keys %{$self->{'ldContainer'}}){
	my %ld_value;  #contains a single ld value in the container {variation_feature variation_feature d_prime r2 snp_distance_count}
	($variation_feature_id_1, $variation_feature_id_2) =  split /-/,$key_ld; #get the variation_features ids
	if (exists $self->{'ldContainer'}->{$key_ld}->{$self->{'_default_population'}}){
	    #add the information to the ld_value hash
	    $ld_value{'variation1'} = $self->{'variationFeatures'}->{$variation_feature_id_1};
	    $ld_value{'variation2'} = $self->{'variationFeatures'}->{$variation_feature_id_2};
	    $ld_value{'d_prime'} = $self->{'ldContainer'}->{$key_ld}->{$self->{'_default_population'}}->{'d_prime'};
	    $ld_value{'r2'} = $self->{'ldContainer'}->{$key_ld}->{$self->{'_default_population'}}->{'r2'};
	    $ld_value{'sample_count'} = $self->{'ldContainer'}->{$key_ld}->{$self->{'_default_population'}}->{'sample_count'};
	    $ld_value{'sample_id'} = $self->{'_default_population'};
	    push @ld_values, \%ld_value;
	}
    }
    return \@ld_values;
}


=head2 get_all_r_square_values

    Example     : $r_square_values = $obj->get_all_r_square_values();
    Description : Get all r_square values contained in the LDFeatureContainer object
    ReturnType  : reference to list of [{variation1=>Bio::EnsEMBL::Variation::VariationFeature, variation2=>Bio::EnsEMBL::Variation::VariationFeature, r2=>r2, sample_id=>population_sample_id}]
    Exceptions  : no exceptions
    Caller      : general
    Status      : At Risk

=cut


sub get_all_r_square_values{
    my $self = shift;
    my @r_squares; #contains ALL the r2 values in the container

    #the keys in the ldContainer hash
    my $variation_feature_id_1;
    my $variation_feature_id_2;

   if (! defined $self->{'_default_population'}){
       $self->{'_default_population'} = $self->_get_major_population;
   }
    foreach my $key_ld (keys %{$self->{'ldContainer'}}){
	my %r_square;  #contains a single r2 value in the container {variation_feature r2 population_id}
	($variation_feature_id_1, $variation_feature_id_2) =  split /-/,$key_ld; #get the variation_features ids
	if (exists $self->{'ldContainer'}->{$key_ld}->{$self->{'_default_population'}}){
	    $r_square{'variation1'} = $self->{'variationFeatures'}->{$variation_feature_id_1};
	    $r_square{'variation2'} = $self->{'variationFeatures'}->{$variation_feature_id_2};
	    $r_square{'r2'} = $self->{'ldContainer'}->{$key_ld}->{$self->{'_default_population'}}->{'r2'};
	    $r_square{'sample_id'} = $self->{'_default_population'};	    
	    #and add all the ld information to the final list
	    push @r_squares, \%r_square;
	}
    }
    return \@r_squares;
}

=head2 get_all_d_prime_values

    Example     : $d_prime_values = $obj->get_all_d_prime_values();
    Description : Get all d_prime values contained in the LDFeatureContainer object
    ReturnType  : reference to list of [{variation1=>Bio::EnsEMBL::Variation::VariationFeature, variation2=>Bio::EnsEMBL::Variation::VariationFeature, d_prime=>d_prime, sample_id=>population_sample_id}]
    Exceptions  : no exceptions
    Caller      : general
    Status      : At Risk

=cut


sub get_all_d_prime_values{
   my $self = shift;
    my @d_primes; #contains ALL the d_prime values in the container

    #the keys in the ldContainer hash
    my $variation_feature_id_1;
    my $variation_feature_id_2;

   if (! defined $self->{'_default_population'}){
       $self->{'_default_population'} = $self->_get_major_population;
   }
   foreach my $key_ld (keys %{$self->{'ldContainer'}}){
       my %d_prime;  #contains a single d_prime value in the container {variation_feature d_prime population_id}
       ($variation_feature_id_1, $variation_feature_id_2) =  split /-/,$key_ld; #get the variation_features ids
       #add the information to the ld_value array
       if (exists $self->{'ldContainer'}->{$key_ld}->{$self->{'_default_population'}}){
	   $d_prime{'variation1'} = $self->{'variationFeatures'}->{$variation_feature_id_1};
	   $d_prime{'variation2'} = $self->{'variationFeatures'}->{$variation_feature_id_2};
	   $d_prime{'d_prime'} = $self->{'ldContainer'}->{$key_ld}->{$self->{'_default_population'}}->{'d_prime'};
	   $d_prime{'sample_id'} = $self->{'_default_population'};
	   #and add all the ld information to the final list if exists the value
	   push @d_primes, \%d_prime;
       }

   }
   return \@d_primes;
}

=head2 get_all_populations

    Arg [1]     : (optional) Bio::EnsEMBL::Variation::VariationFeature $variationFeature
    Arg [2]     : (optional) Bio::EnsEMBL::Variation::VariationFeature $variationFeature
    Example     : $populations = $obj->get_all_populations($vf1,$vf2);
    Description : If no arguments provided, returns ALL the populations present in the container. When 2 variation features provided, returns the 
    population/populations where these variation features occurs
    ReturnType  : ref to an array of int
    Exceptions  : throw on incorrect arguments
    Caller      : general
    Status      : At Risk

=cut

sub get_all_populations{
    my $self = shift;
    my $variation_feature_1 = shift;
    my $variation_feature_2 = shift;
    my %populations;
    my @populations;    
    my $key;

    #if no variation provided, return ALL the populations in the container
    if (! defined($variation_feature_1) && ! defined($variation_feature_2)){
	foreach my $key (keys %{$self->{'ldContainer'}}){
	    map {$populations{$_}++} keys %{$self->{'ldContainer'}->{$key}};
	}
	@populations = keys %populations;
    }
    else{
	#first, check if both arguments have been properly provided
	if (!ref($variation_feature_1) || !$variation_feature_1->isa('Bio::EnsEMBL::Variation::VariationFeature') || !ref($variation_feature_2) || !$variation_feature_2->isa('Bio::EnsEMBL::Variation::VariationFeature')){
	    throw("Bio::EnsEMBL::Variation::VariationFeature arguments expected");
	}
	#if the variation_features are correct, return the list of populations
	else{
	    #find out the key in the ld Hash: vf1 - vf2 or vf2 - vf1
	    if (exists $self->{'ldContainer'}->{$variation_feature_1->dbID() . '-' . $variation_feature_2->dbID()}){
		$key = $variation_feature_1->dbID() . '-' . $variation_feature_2->dbID();
	    }
	    else{
		$key = $variation_feature_2->dbID() . '-' . $variation_feature_1->dbID();
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
1;

