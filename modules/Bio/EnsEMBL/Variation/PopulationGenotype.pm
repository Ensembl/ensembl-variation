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
  Status     : At Risk

=cut

sub new {
  my $class = shift;

  my ($dbID, $adaptor, $allele1, $allele2, $var, $pop, $freq, $count, $ss_id) =
    rearrange([qw(dbID adaptor allele1 allele2 
                  variation population frequency count subsnp)],@_);

  if(defined($var) &&
     (!ref($var) || !$var->isa('Bio::EnsEMBL::Variation::Variation'))) {
    throw("Bio::EnsEMBL::Variation::Variation argument expected");
  }

  if(defined($pop) &&
     (!ref($pop) || !$pop->isa('Bio::EnsEMBL::Variation::Population'))) {
    throw("Bio::EnsEMBL::Variation::Population argument expected");
  }
  
  # set subsnp_id to undefined if it's 0 in the DB
  $ss_id = undef if defined($ss_id) && $ss_id == 0;
  
  # add ss to the subsnp_id
  $ss_id = 'ss'.$ss_id if defined $ss_id && $ss_id !~ /^ss/;

  return bless {'dbID'    => $dbID,
                'adaptor' => $adaptor,
                'allele1' => $allele1,
                'allele2' => $allele2,
                'variation' => $var,
                'population' => $pop,
                'frequency' => $freq,
                'count'   => $count,
                'subsnp' => $ss_id}, $class;
}




=head2 population

  Arg [1]    : (optional) Bio::EnsEMBL::Variation::Population $pop
  Example    : $pop = $pop_genotype->population();
  Description: Getter/Setter for the population associated with this genotype
  Returntype : Bio::EnsEMBL::Variation::Population
  Exceptions : throw on bad argument
  Caller     : general
  Status     : At Risk

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
  Status     : At Risk

=cut

sub frequency{
  my $self = shift;
  return $self->{'frequency'} = shift if(@_);
  return $self->{'frequency'};
}

=head2 count

  Arg [1]    : int $count (optional)
               The new value to set the count attribute to
  Example    : $frequency = $pop_gtype->count()
  Description: Getter/Setter for the observed count of this genotype
               within its associated population.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub count{
  my $self = shift;
  return $self->{'count'} = shift if(@_);
  return $self->{'count'};
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
    my $v = shift;
    if(defined($v) &&
       (!ref($v) || !$v->isa('Bio::EnsEMBL::Variation::Variation'))) {
      throw('Bio::EnsEMBL::Variation::Variation argument expected.');
    }
    return $self->{'variation'} = $v;
  }
  return $self->{'variation'};
}



=head2 subsnp

  Arg [1]    : string $newval (optional) 
               The new value to set the subsnp attribute to
  Example    : print $a->subsnp();
  Description: Getter/Setter for the subsnp attribute.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub subsnp{
  my $self = shift;
  return $self->{'subsnp'} = shift if(@_);
  return $self->{'subsnp'};
}



=head2 subsnp_handle

  Arg [1]    : string $newval (optional) 
               The new value to set the subsnp_handle attribute to
  Example    : print $a->subsnp_handle();
  Description: Getter/Setter for the subsnp_handle attribute.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub subsnp_handle{
  my $self = shift;
  
  # if changing handle
  if(@_) {
    return $self->{'subsnp_handle'} = shift;
  }
  
  # if not already defined, retrieve from the database
  if(!defined $self->{'subsnp_handle'}) {
    
    # check if the subsnp is useable and the db exists
    if(defined ($self->{'subsnp'}) && defined ($self->{'adaptor'})) {
      my $ss = $self->subsnp();
      
      # get rid of the ss from the beginning
      $ss =~ s/^ss//g;
      
      my $sth = $self->{'adaptor'}->dbc->prepare(qq/SELECT handle FROM subsnp_handle WHERE subsnp_id = ?;/);
      $sth->execute($ss);
      
      my $handle = $sth->fetchrow_arrayref->[0];
      
      return $self->{'subsnp_handle'} = $handle;
    }
  }
  
  return $self->{'subsnp_handle'};
}

1;
