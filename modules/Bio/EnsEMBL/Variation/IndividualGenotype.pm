# Ensembl module for Bio::EnsEMBL::Variation::IndividualGenotype
#
# Copyright (c) 2004 Ensembl
#


=head1 NAME

Bio::EnsEMBL::Variation::IndividualGenotype - Module representing the genotype
of a single individual at a single locus

=head1 SYNOPSIS

    print $genotype->variation()->name(), "\n";
    print $genotype->allele1(), '/', $genotype->allele2(), "\n";
    print $genotype->individual()->name(), "\n";

=head1 DESCRIPTION

This is a class representing the genotype of a single diploid individual at
a specific locus.

=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::Variation::IndividualGenotype;

use Bio::EnsEMBL::Variation::Genotype;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw deprecate warning);

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Variation::Genotype);



=head2 new

  Arg [-adaptor] :
    Bio::EnsEMBL::Variation::DBSQL::IndividualAdaptor
  Arg [-allele1] :
    string - One of the two alleles defining this genotype
  Arg [-allele2] :
    string - One of the two alleles defining this genotype
  Arg [-variation] :
    Bio::EnsEMBL::Variation::Variation - The variation associated with this
    genotype
  Arg [-individual] :
    Bio::EnsEMBL::Individual - The individual this genotype is for.
  Example    : $ind_genotype = Bio:EnsEMBL::Variation::IndividualGenotype->new
                   (-allele1 => 'A',
                    -allele2 => 'T',
                    -variation => $variation,
                    -individual => $ind);
  Description: Constructor.  Instantiates an IndividualGenotype object.
  Returntype : Bio::EnsEMBL::Variation::IndividualGenotype
  Exceptions : throw on bad argument
  Caller     : general

=cut

sub new {
  my $class = shift;

  my ($adaptor, $allele1, $allele2, $var, $ind) =
    rearrange([qw(adaptor allele1 allele2 variation individual)],@_);

  if(defined($var) &&
     (!ref($var) || !$var->isa('Bio::EnsEMBL::Variation::Variation'))) {
    throw("Bio::EnsEMBL::Variation::Variation argument expected");
  }

  if(defined($ind) &&
     (!ref($ind) || !$ind->isa('Bio::EnsEMBL::Variation::Individual'))) {
    throw("Bio::EnsEMBL::Variation::Individual argument expected");
  }

  return bless {'adaptor' => $adaptor,
                'allele1' => $allele1,
                'allele2' => $allele2,
                'variation' => $var,
                'individual' => $ind}, $class;
}




=head2 individual

  Arg [1]    : (optional) Bio::EnsEMBL::Variation::Individual $ind
  Example    : $ind = $ind_genotype->individual();
  Description: Getter/Setter for the individual associated with this genotype
  Returntype : Bio::EnsEMBL::Variation::Individual
  Exceptions : throw on bad argument
  Caller     : general

=cut


sub individual {
  my $self = shift;
  if(@_) {
    my $ind = shift;
    if(defined($ind) &&
       (!ref($ind) || !$ind->isa('Bio::EnsEMBL::Variation::Individual'))) {
      throw('Bio::EnsEMBL::Variation::Individual argument expected');
    }
    return $self->{'individual'} = $ind;
  }
  return $self->{'individual'};
}

1;
