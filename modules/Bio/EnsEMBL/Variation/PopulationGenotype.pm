# Ensembl module for Bio::EnsEMBL::Variation::PopulationGenotype
#
# Copyright (c) 2004 Ensembl
#


=head1 NAME

Bio::EnsEMBL::Variation::PopulationGenotype - Module for a genotype
represented in a population.

=head1 SYNOPSIS

    print $genotype->variation()->name(), "\n";
    print $genotype->allele1(), '/', $genotype->allele2(), "\n";
    print $genotype->frequency(), "\n";
    print $genotype->population()->name(), "\n";

=head1 DESCRIPTION

This class represents a genotype which is present in a population.

=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::Variation::PopulationGenotype;

use Bio::EnsEMBL::Variation::Genotype;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw deprecate warning);

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Variation::Genotype);



=head2 new

  Arg [-dbID] :
    int - unique internal identifier
  Arg [-adaptor] :
    Bio::EnsEMBL::Variation::DBSQL::PopulationAdaptor
  Arg [-allele1] :
    string - One of the two alleles defining this genotype
  Arg [-allele2] :
    string - One of the two alleles defining this genotype
  Arg [-variation] :
    Bio::EnsEMBL::Variation::Variation - The variation associated with this
    genotype
  Arg [-population] :
    Bio::EnsEMBL::Population - The population this genotype is for.
  Arg [-frequency] :
    int - the frequency this genotype occurs in this population
  Example    : $pop_genotype = Bio:EnsEMBL::Variation::PopulationGenotype->new
                   (-allele1 => 'A',
                    -allele2 => 'T',
                    -variation => $variation,
                    -population => $pop
                    -frequency  => 0.87);
  Description: Constructor.  Instantiates a PopulationGenotype object.
  Returntype : Bio::EnsEMBL::Variation::PopulationGenotype
  Exceptions : throw on bad argument
  Caller     : general

=cut

sub new {
  my $class = shift;

  my ($dbID, $adaptor, $allele1, $allele2, $var, $pop, $freq) =
    rearrange([qw(dbID adaptor allele1 allele2 
                  variation population frequency)],@_);

  if(defined($var) &&
     (!ref($var) || !$var->isa('Bio::EnsEMBL::Variation::Variation'))) {
    throw("Bio::EnsEMBL::Variation::Variation argument expected");
  }

  if(defined($pop) &&
     (!ref($pop) || !$pop->isa('Bio::EnsEMBL::Variation::Population'))) {
    throw("Bio::EnsEMBL::Variation::Population argument expected");
  }

  return bless {'dbID'    => $dbID,
                'adaptor' => $adaptor,
                'allele1' => $allele1,
                'allele2' => $allele2,
                'variation' => $var,
                'population' => $pop,
                'frequency' => $freq}, $class;
}




=head2 population

  Arg [1]    : (optional) Bio::EnsEMBL::Variation::Population $pop
  Example    : $pop = $pop_genotype->population();
  Description: Getter/Setter for the population associated with this genotype
  Returntype : Bio::EnsEMBL::Variation::Population
  Exceptions : throw on bad argument
  Caller     : general

=cut


sub population {
  my $self = shift;
  if(@_) {
    my $pop = shift;
    if(defined($pop) &&
       (!ref($pop) || !$pop->isa('Bio::EnsEMBL::Variation::Population'))) {
      throw('Bio::EnsEMBL::Variation::Population argument expected');
    }
    return $self->{'population'} = $pop;
  }
  return $self->{'population'};
}




=head2 frequency

  Arg [1]    : string $freq (optional)
               The new value to set the frequency attribute to
  Example    : $frequency = $pop_gtype->frequency()
  Description: Getter/Setter for the frequency of occurance of this genotype
               within its associated population.
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub frequency{
  my $self = shift;
  return $self->{'frequency'} = shift if(@_);
  return $self->{'frequency'};
}

1;
