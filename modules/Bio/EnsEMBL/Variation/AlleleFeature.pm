# Ensembl module for Bio::EnsEMBL::Variation::AlleleFeature
#
# Copyright (c) 2005 Ensembl
#


=head1 NAME

Bio::EnsEMBL::Variation::AlleleFeature - A genomic position for an allele in a population.

=head1 SYNOPSIS

    # Allele feature representing a single nucleotide polymorphism
    $af = Bio::EnsEMBL::Variation::AlleleFeature->new
       (-start   => 100,
        -end     => 100,
        -strand  => 1,
        -slice   => $slice,
        -allele_string => 'A',
        -variation_name => 'rs635421',
        -variation => $v);
    ...

    # a allele feature is like any other ensembl feature, can be
    # transformed etc.
    $af = $af->transform('supercontig');

    print $af->start(), "-", $af->end(), '(', $af->strand(), ')', "\n";

    print $af->name(), ":", $af->allele_string();

    # Get the Variation object which this feature represents the genomic
    # position of. If not already retrieved from the DB, this will be
    # transparently lazy-loaded
    my $v = $af->variation();

=head1 DESCRIPTION

This is a class representing the genomic position of a allele in a population
from the ensembl-variation database.  The actual variation information is
represented by an associated Bio::EnsEMBL::Variation::Variation object. Some
of the information has been denormalized and is available on the feature for
speed purposes.  A AlleleFeature behaves as any other Ensembl feature.
See B<Bio::EnsEMBL::Feature> and B<Bio::EnsEMBL::Variation::Variation>.

=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::AlleleFeature;

use Bio::EnsEMBL::Feature;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument  qw(rearrange);


our @ISA = ('Bio::EnsEMBL::Feature');

=head2 new

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

  Arg [-VARIATION] :
    int - the variation object associated with this feature.

  Arg [-VARIATION_ID] :
    int - the internal id of the variation object associated with this
    identifier. This may be provided instead of a variation object so that
    the variation may be lazy-loaded from the database on demand.
    
  Arg [-POPULATION_ID] :
    int - the internal id of the population object associated with this
    identifier. This may be provided instead of the population object so that
    the population may be lazy-loaded from the database on demand.

  Arg [-ALLELE_STRING] :
    string - the allele for this AlleleFeature object.

  Example    :
    $af = Bio::EnsEMBL::Variation::AlleleFeature->new
       (-start   => 100,
        -end     => 100,
        -strand  => 1,
        -slice   => $slice,
        -allele_string => 'A',
        -variation_name => 'rs635421',
	-population  => $p,
        -variation => $v);

  Description: Constructor. Instantiates a new AlleleFeature object.
  Returntype : Bio::EnsEMBL::Variation::AlleleFeature
  Exceptions : none
  Caller     : general

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);
  my ($allele, $var_name, $variation, $variation_id,$population, $population_id) =
    rearrange([qw(ALLELE_STRING VARIATION_NAME 
                  VARIATION VARIATION_ID POPULATION POPULATION_ID)], @_);

  $self->{'allele'}           = $allele;
  $self->{'variation_name'}   = $var_name;
  $self->{'variation'}        = $variation;
  $self->{'_variation_id'}    = $variation_id;
  $self->{'population'}       = $population;
  $self->{'_population_id'}   = $population_id;

  return $self;
}



sub new_fast {
  my $class = shift;
  my $hashref = shift;
  return bless $hashref, $class;
}


=head2 allele_string

  Arg [1]    : string $newval (optional)
               The new value to set the allele attribute to
  Example    : $allele = $obj->allele_string()
  Description: Getter/Setter for the allele attribute.
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub allele_string{
  my $self = shift;
  return $self->{'allele_string'} = shift if(@_);
  return $self->{'allele_string'};
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

=head2 variation

  Arg [1]    : (optional) Bio::EnsEMBL::Variation::Variation $variation
  Example    : $v = $af->variation();
  Description: Getter/Setter for the variation associated with this feature.
               If not set, and this AlleleFeature has an associated adaptor
               an attempt will be made to lazy-load the variation from the
               database.
  Returntype : Bio::EnsEMBL::Variation::Variation
  Exceptions : throw on incorrect argument
  Caller     : general

=cut

sub variation {
  my $self = shift;

  if(@_) {
    if(!ref($_[0]) || !$_[0]->isa('Bio::EnsEMBL::Variation::Variation')) {
      throw("Bio::EnsEMBL::Variation::Variation argument expected");
    }
    $self->{'variation'} = shift;
  }
  elsif(!defined($self->{'variation'}) && $self->{'adaptor'} &&
        defined($self->{'_variation_id'})) {
    # lazy-load from database on demand
    my $va = $self->{'adaptor'}->db()->get_VariationAdaptor();
    $self->{'variation'} = $va->fetch_by_dbID($self->{'_variation_id'});
  }

  return $self->{'variation'};
}

=head2 population

  Arg [1]    : (optional) Bio::EnsEMBL::Variation::Population $population
  Example    : $p = $af->population();
  Description: Getter/Setter for the population associated with this feature.
               If not set, and this AlleleFeature has an associated adaptor
               an attempt will be made to lazy-load the population from the
               database.
  Returntype : Bio::EnsEMBL::Variation::Population
  Exceptions : throw on incorrect argument
  Caller     : general

=cut

sub population {
  my $self = shift;

  if(@_) {
    if(!ref($_[0]) || !$_[0]->isa('Bio::EnsEMBL::Variation::Population')) {
      throw("Bio::EnsEMBL::Variation::Population argument expected");
    }
    $self->{'population'} = shift;
  }
  elsif(!defined($self->{'population'}) && $self->{'adaptor'} &&
        defined($self->{'_population_id'})) {
    # lazy-load from database on demand
    my $pa = $self->{'adaptor'}->db()->get_PopulationAdaptor();
    $self->{'population'} = $pa->fetch_by_dbID($self->{'_population_id'});
  }

  return $self->{'population'};
}

=head2 apply_edit
    
    Arg [1]    : reference to string $seqref
    Arg [2]    : int $start of the seq_ref
    Example    : $sequence = 'ACTGAATATTTAAGGCA';
               $af->apply_edit(\$sequence,$start);
               print $sequence, "\n";
    Description: Applies this AlleleFeature directly to a sequence which is
               passed by reference.  
               If either the start or end of this AlleleFeature are not defined
               this function will not do anything to the passed sequence.
    Returntype : reference to the same sequence that was passed in
    Exceptions : none
    Caller     : Slice

=cut

sub apply_edit  {

  my $self   = shift;
  my $seqref = shift;

  if(ref($seqref) ne 'SCALAR') {
    throw("Reference to scalar argument expected");
  }

  if(!defined($self->{'start'}) || !defined($self->{'end'})) {
    return $seqref;
  }


  my $len = $self->length;
  substr($$seqref, $self->{'start'}-1, $len) = $self->{'allele_string'} if ($self->{'allele_string'} ne '-'); 
  substr($$seqref, $self->{'start'}-1, $len) = '' if ($self->{'allele_string'} eq '-');
  return $seqref;

}

=head2 length_diff

  Arg [1]    : none
  Example    : my $diff = $af->length_diff();
  Description: Returns the difference in length caused by applying this
               AlleleFeature to a sequence.  This may be be negative (deletion),
               positive (insertion) or 0 (replacement).
               If either start or end are not defined 0 is returned.
  Returntype : int
  Exceptions : none
  Caller     : general

=cut

sub length_diff  {

  my $self = shift;

  return 0 if(!defined($self->{'end'}) || !defined($self->{'start'}));

  return length($self->{'allele_string'}) - ($self->{'end'} - $self->{'start'} + 1) if ($self->{'allele_string'} ne '-'); 
  return 0 - ($self->{'end'} - $self->{'start'} + 1) if ($self->{'allele_string'} eq '-'); 

}

1;
