# Ensembl module for Bio::EnsEMBL::Variation::VariationFeature
#
# Copyright (c) 2004 Ensembl
#


=head1 NAME

Bio::EnsEMBL::Variation::VariationFeature - A genomic position for a nucleotide variation.

=head1 SYNOPSIS

    # Variation feature representing a single nucleotide polymorphism
    $vf = Bio::EnsEMBL::Variation::VariationFeature->new
       (-start   => 100,
        -end     => 100,
        -strand  => 1,
        -slice   => $slice,
        -allele_string => 'A/T',
        -variation_name => 'rs635421',
        -map_weight  => 1,
        -variation => $v);

    # Variation feature representing a 2bp insertion
    $vf = Bio::EnsEMBL::Variation::VariationFeature->new
       (-start   => 1522,
        -end     => 1521, # end = start-1 for insert
        -strand  => -1,
        -slice   => $slice,
        -allele_string => '-/AA',
        -variation_name => 'rs12111',
        -map_weight  => 1,
        -variation => $v2);

    ...

    # a variation feature is like any other ensembl feature, can be
    # transformed etc.
    $vf = $vf->transform('supercontig');

    print $vf->start(), "-", $vf->end(), '(', $vf->strand(), ')', "\n";

    print $vf->name(), ":", $vf->allele_string();

    # Get the Variation object which this feature represents the genomic
    # position of. If not already retrieved from the DB, this will be
    # transparently lazy-loaded
    my $v = $vf->get_Variation();

=head1 DESCRIPTION

This is a class representing the genomic position of a nucleotide variation
from the ensembl-variation database.  The actual variation information is
represented by an associated Bio::EnsEMBL::Variation::Variation object. Some
of the information has been denormalized and is available on the feature for
speed purposes.  A VariationFeature behaves as any other Ensembl feature.
See B<Bio::EnsEMBL::Feature> and B<Bio::EnsEMBL::Variation::Variation>.

=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::Variation::VariationFeature;


=head2 new

  Arg [-dbID] :
    see superclass constructor

  Arg [-ADAPTOR] :
    see superclass constructor

  Arg [-ADAPTOR] :
    see superclass constructor

  Arg [-START] :
    see superclass constructor
  Arg [-END] :
    see superclass constructor

  Arg [-STRAND] :
    see superclass constructor

  Arg [-SLICE] :
    see superclass constructor

  Arg [-VARIATION_NAME] :
    string - the name of the variation this feature is for (denormalisation
    from Variation object).

  Arg [-MAP_WEIGHT] :
    int - the number of times that the variation associated with this feature
    has hit the genome. If this was the only feature associated with this
    variation_feature the map_weight would be 1.

  Example    :
    $vf = Bio::EnsEMBL::Variation::VariationFeature->new
       (-start   => 100,
        -end     => 100,
        -strand  => 1,
        -slice   => $slice,
        -allele_string => 'A/T',
        -variation_name => 'rs635421',
        -map_weight  => 1,
        -variation => $v);

  Description: Constructor. Instantiates a new VariationFeature object.
  Returntype : Bio::EnsEMBL::Variation::Variation
  Exceptions : none
  Caller     : general

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);

  my ($allele_str, $var_name, $map_weight, $variation) =
    rearrange([qw(ALLELE_STRING VARIATION_NAME MAP_WEIGHT VARIATION)], @_);

  $self->{'allele_string'}  = $allele_str;
  $self->{'variation_name'} = $var_name;
  $self->{'map_weight'}     = $map_weight;
  $self->{'variation'}      = $variation;

  return $self;
}



=head2 allele_string

  Arg [1]    : string $newval (optional)
               The new value to set the allele_string attribute to
  Example    : $allele_string = $obj->allele_string()
  Description: Getter/Setter for the allele_string attribute.
               The allele_string is a '/' demimited string representing the
               alleles associated with this features variation.
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub allele_string{
  my $self = shift;
  return $self->{'allele_string'} = shift if(@_);
  return $self->{'allele_string'};
}



=head2 display_id

  Arg [1]    : none
  Example    : print $vf->display_id(), "\n";
  Description: Returns the 'display' identifier for this feature. For
               VariationFeatures this is simply the name of the variation
               it is associated with.
  Returntype : string
  Exceptions : none
  Caller     : webcode

=cut

sub display_id {
  my $self = shift;
  return $self->{'variation_name'} || '';
}



=head2 variation_name

  Arg [1]    : string $newval (optional)
               The new value to set the variation_name attribute to
  Example    : $variation_name = $obj->variation_name()
  Description: Getter/Setter for the variation_name attribute.  This is the
               name of the variation associated with this feature.
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub variation_name{
  my $self = shift;
  return $self->{'variation_name'} = shift if(@_);
  return $self->{'variation_name'};
}



=head2 map_weight

  Arg [1]    : int $newval (optional) 
               The new value to set the map_weight attribute to
  Example    : $map_weight = $obj->map_weight()
  Description: Getter/Setter for the map_weight attribute. The map_weight
               is the number of times this features variation was mapped to
               the genome.
  Returntype : int
  Exceptions : none
  Caller     : general

=cut

sub map_weight{
  my $self = shift;
  return $self->{'map_weight'} = shift if(@_);
  return $self->{'map_weight'};
}



=head2 get_Variation

  Arg [1]    : none
  Example    : $v = $vf->get_Variation();
  Description: Retrieves the variation associated with this feature.
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub get_Variation {
  my $self = shift;

  ### TODO: this should be lazy-loaded on demand
  return $self->{'variation'};
}


1;
