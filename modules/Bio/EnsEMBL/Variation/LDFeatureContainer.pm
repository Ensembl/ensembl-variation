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
        -LdContainer => $ldhash);

    ...

    #get the d_prime values for a certain variation_features
    D_prime = get_d_prime($variation_feature_id_1,$variation_feature_id_2);
    #get the list of variation in the container
    $variations = get_variations();

    ...

=head1 DESCRIPTION

This is a class representing the LD information for a certain region
from the ensembl-variation database.
See B<Bio::EnsEMBL::Variation::Variation>.

=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::LDFeatureContainer;

use Bio::EnsEMBL::Storable;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument  qw(rearrange);

use vars qw(@ISA);
use Data::Dumper;

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
    $ldContainer = Bio::EnsEMBL::Variation::LDFeatureContainer->new
       (-name => 'chr1'
        -ldContainer => {'variation_feature_1-variation_feature_2' => { 'population_id_1' =>
	                                                                          { 'Dprime' => 0.5,
										    'r2'      => 0.421,
										    'snp_distance_count' => 5,
										    'sample_count' => 120
										    },
									'population_id_2' => 
 								                  { 'Dprime' => 0.3,
										    'r2'     => 0.321,
										    'snp_distance_count' => 3,
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
    Arg [3]     : int - population you want to get the r_square value
    Example     : $r_square = $obj->get_r_square($vf1,$vf2,$population_id);
    Description : Get the r_square value for a pair of variation features in the given population. If no population is provided,
    return the r_square for the population with more sample counts (in case more than 1)
    ReturnType  : float
    Exceptions  : throw on incorrect arguments
    Caller      : general
=cut

sub get_r_square{
    my $self = shift;
    my $variation_feature_1 = shift;
    my $variation_feature_2 = shift;
    my $population_id = shift;
    
    $population_id ||= 0;
    my $r_square;
    my $key;
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
	    if (exists $self->{'ldContainer'}->{$key}->{$population_id}){
		$r_square = $self->{'ldContainer'}->{$key}->{$population_id}->{'r2'}
	    }
	    else{
		#there was no population provided, return the r_square where the population has a bigger sample_count
		$population_id = $self->_get_major_population($self->{'ldContainer'}->{$key});
		$r_square = $self->{'ldContainer'}->{$key}->{$population_id}->{'r2'};
	    }
	}
	    
    }
    return $r_square;
}

=head2 get_d_prime

    Arg [1]     : Bio::EnsEMBL::Variation::VariationFeature $variationFeature
    Arg [2]     : Bio::EnsEMBL::Variation::VariationFeature $variationFeature
    Arg [3]     : int - population you want to get the d_prime value
    Example     : $d_prime = $obj->get_d_prime($vf1,$vf2,$population_id);
    Description : Get the d_prime value for a pair of variation features for a known or unknown population
    ReturnType  : float
    Exceptions  : throw on incorrect arguments
    Caller      : general
=cut

sub get_d_prime{
    my $self = shift;
    my $variation_feature_1 = shift;
    my $variation_feature_2 = shift;
    my $population_id = shift;
    
    $population_id ||= 0;
    my $d_prime;
    my $key;
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
	    if (exists $self->{'ldContainer'}->{$key}->{$population_id}){
		$d_prime = $self->{'ldContainer'}->{$key}->{$population_id}->{'Dprime'}
	    }
	    else{
		#there was no population provided, return the r_square where the population has a bigger sample_count
		$population_id = $self->_get_major_population($self->{'ldContainer'}->{$key});
		$d_prime = $self->{'ldContainer'}->{$key}->{$population_id}->{'Dprime'};
	    }
	}
	    
    }
    return $d_prime;
}


=head2 get_all_ld_values

    Example     : $ld_values = $obj->get_all_ld_values();
    Description : Get all the information contained in the LDFeatureContainer object
    ReturnType  : reference to list of [Bio::EnsEMBL::Variation::VariationFeature Bio::EnsEMBL::Variation::VariationFeature Dprime r2 snp_distance_count sample_count population_id]
    Exceptions  : thrown when there is not valid db connection
    Caller      : general
=cut


sub get_all_ld_values{
    my $self = shift;
    my @ld_values; #contains ALL the ld values in the container

    #the keys in the ldContainer hash
    my $variation_feature_id_1;
    my $variation_feature_id_2;

    my $population_id;

    foreach my $key_ld (keys %{$self->{'ldContainer'}}){
	my @ld_value;  #contains a single ld value in the container [variation_feature variation_feature Dprime r2 snp_distance_count]
	($variation_feature_id_1, $variation_feature_id_2) =  split /-/,$key_ld; #get the variation_features ids
	#add the information to the ld_value array
	push @ld_value, $self->{'variationFeatures'}->{$variation_feature_id_1}, $self->{'variationFeatures'}->{$variation_feature_id_2};
	#and add the ld information: dprime, r2, snp_distance_count and sample_count for the major population (population where the sample count is higher
	$population_id = $self->_get_major_population($self->{'ldContainer'}->{$key_ld});

	push @ld_value, $self->{'ldContainer'}->{$key_ld}->{$population_id}->{'Dprime'}, $self->{'ldContainer'}->{$key_ld}->{$population_id}->{'r2'},$self->{'ldContainer'}->{$key_ld}->{$population_id}->{'snp_distance_count'}, $self->{'ldContainer'}->{$key_ld}->{$population_id}->{'sample_count'}, $population_id;
	#and add all the ld information to the final list
	push @ld_values, \@ld_value;
    }
    return \@ld_values;
}


=head2 get_all_r_square_values

    Example     : $r_square_values = $obj->get_all_r_square_values();
    Description : Get all r_square values contained in the LDFeatureContainer object
    ReturnType  : reference to list of [Bio::EnsEMBL::Variation::VariationFeature Bio::EnsEMBL::Variation::VariationFeature r2 population_id]
    Exceptions  : thrown when there is not valid db connection
    Caller      : general
=cut


sub get_all_r_square_values{
    my $self = shift;
    my @r_squares; #contains ALL the r2 values in the container

    #the keys in the ldContainer hash
    my $variation_feature_id_1;
    my $variation_feature_id_2;

    my $population_id;

    foreach my $key_ld (keys %{$self->{'ldContainer'}}){
	my @r_square;  #contains a single r2 value in the container [variation_feature r2 population_id]
	($variation_feature_id_1, $variation_feature_id_2) =  split /-/,$key_ld; #get the variation_features ids
	#add the information to the ld_value array
	push @r_square, $self->{'variationFeatures'}->{$variation_feature_id_1}, $self->{'variationFeatures'}->{$variation_feature_id_2};
	#and add the r2 information
	$population_id = $self->_get_major_population($self->{'ldContainer'}->{$key_ld});

	push @r_square, $self->{'ldContainer'}->{$key_ld}->{$population_id}->{'r2'}, $population_id;
	#and add all the ld information to the final list
	push @r_squares, \@r_square;
    }
    return \@r_squares;
}

=head2 get_all_d_prime_values

    Example     : $d_prime_values = $obj->get_all_d_prime_values();
    Description : Get all d_prime values contained in the LDFeatureContainer object
    ReturnType  : reference to list of [Bio::EnsEMBL::Variation::VariationFeature Bio::EnsEMBL::Variation::VariationFeature d_prime population_id]
    Exceptions  : thrown when there is not valid db connection
    Caller      : general
=cut


sub get_all_d_prime_values{
   my $self = shift;
    my @d_primes; #contains ALL the d_prime values in the container

    #the keys in the ldContainer hash
    my $variation_feature_id_1;
    my $variation_feature_id_2;

    my $population_id;

    foreach my $key_ld (keys %{$self->{'ldContainer'}}){
	my @d_prime;  #contains a single r2 value in the container [variation_feature r2 population_id]
	($variation_feature_id_1, $variation_feature_id_2) =  split /-/,$key_ld; #get the variation_features ids
	#add the information to the ld_value array
	push @d_prime, $self->{'variationFeatures'}->{$variation_feature_id_1}, $self->{'variationFeatures'}->{$variation_feature_id_2};
	#and add the r2 information
	$population_id = $self->_get_major_population($self->{'ldContainer'}->{$key_ld});

	push @d_prime, $self->{'ldContainer'}->{$key_ld}->{$population_id}->{'Dprime'}, $population_id;
	#and add all the ld information to the final list
	push @d_primes, \@d_prime;
    }
    return \@d_primes;
}


sub _get_major_population{
    my $self = shift;
    my $ld_population = shift;
    my $population_id;
    my $max_sample_count = 0;
    foreach my $population (keys %{$ld_population}){
	if ($ld_population->{$population}->{'sample_count'} > $max_sample_count){
	    $max_sample_count = $ld_population->{$population}->{'sample_count'};
	    $population_id = $population;
	}
    }
    return $population_id;
}
1;

