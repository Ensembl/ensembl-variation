=head1 LICENSE

 Copyright (c) 1999-2011 The European Bioinformatics Institute and
 Genome Research Limited.  All rights reserved.

 This software is distributed under a modified Apache license.
 For license details, please see

   http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <dev@ensembl.org>.

 Questions may also be sent to the Ensembl help desk at
 <helpdesk@ensembl.org>.

=cut

# Ensembl module for Bio::EnsEMBL::Variation::IndividualGenotype
#
# Copyright (c) 2004 Ensembl
#


=head1 NAME

Bio::EnsEMBL::Variation::IndividualGenotype- Module representing the genotype
of a single individual at a single position

=head1 SYNOPSIS

    print $genotype->variation()->name(), "\n";
    print $genotype->allele1(), '/', $genotype->allele2(), "\n";
    print $genotype->individual()->name(), "\n";

=head1 DESCRIPTION

This is a class representing the genotype of a single diploid individual at
a specific position

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::Variation::IndividualGenotypeFeature;

use Bio::EnsEMBL::Variation::Genotype;
use Bio::EnsEMBL::Feature;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw deprecate warning);

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Variation::Genotype Bio::EnsEMBL::Feature);



=head2 new

  Arg [-adaptor] :
    Bio::EnsEMBL::Variation::DBSQL::IndividualAdaptor
  Arg [-START] :
    see superclass constructor
  Arg [-END] :
    see superclass constructor
  Arg [-STRAND] :
    see superclass constructor
  Arg [-SLICE] :
    see superclass constructor
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
                   (-start   => 100,
		    -end     => 100,
		    -strand  => 1,
		    -slice   => $slice,
		    -allele1 => 'A',
                    -allele2 => 'T',
                    -variation => $variation,
                    -individual => $ind);
  Description: Constructor.  Instantiates an IndividualGenotype object.
  Returntype : Bio::EnsEMBL::Variation::IndividualGenotype
  Exceptions : throw on bad argument
  Caller     : general
  Status     : At Risk

=cut

sub new {
    my $caller = shift;
    my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);

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

    $self->{'adaptor'} = $adaptor;
    $self->{'allele1'} = $allele1;
    $self->{'allele2'} = $allele2;
    $self->{'individual'} = $ind;
    $self->{'variation'} = $var;
    
    return $self;

}



sub new_fast {
  my $class = shift;
  my $hashref = shift;
  return bless $hashref, $class;
}


=head2 individual

  Arg [1]    : (optional) Bio::EnsEMBL::Variation::Individual $ind
  Example    : $ind = $ind_genotype->individual();
  Description: Getter/Setter for the individual associated with this genotype
  Returntype : Bio::EnsEMBL::Variation::Individual
  Exceptions : throw on bad argument
  Caller     : general
  Status     : At Risk

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

=head2 variation

  Arg [1]    : (optional) Bio::EnsEMBL::Variation::Variation $var
  Example    : $var = $genotype->variation();
  Description: Getter/Setter for the Variation as
  Returntype : Bio::EnsEMBL::Variation::Variation
  Exceptions : throw on bad argument
  Caller     : general
  Status     : At Risk

=cut

sub variation {
  my $self = shift;

  if(@_) {
      #Setter: check wether it is a variation and return it
      my $v = shift;
      if(defined($v) &&
	 (!ref($v) || !$v->isa('Bio::EnsEMBL::Variation::Variation'))) {
	  throw('Bio::EnsEMBL::Variation::Variation argument expected.');
      }
      $self->{'variation'} = $v;
  }
  else{
      if(!defined($self->{'variation'}) && $self->{'adaptor'})    {
	  #lazy-load from database on demand
	  my $vfa = $self->{'adaptor'}->db()->get_VariationFeatureAdaptor();
	  $self->{'variation'} = (shift @{$vfa->fetch_all_by_Slice($self->feature_Slice())})->variation;
      }
  }
  return $self->{'variation'};
}

1;
