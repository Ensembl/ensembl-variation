# Ensembl module for Bio::EnsEMBL::Variation::Allele
#
# Copyright (c) 2004 Ensembl
#


=head1 NAME

Bio::EnsEMBL::Variation::Allele - A single allele of a nucleotide variation.

=head1 SYNOPSIS

    $allele = Bio::EnsEMBL::Variation::Allele->new
       (-allele => 'A',
        -frequency => 0.85,
        -population => $population);

    $delete = Bio::EnsEMBL::Variation::Allele->new
       (-allele => '-',
        -frequency => 0.15,
        -population => $population);

    ...

    $astr = $a->allele();
    $pop  = $a->population();
    $freq = $a->frequency();

    print $a->allele();
    if($a->populaton) {
       print " found in population ", $allele->population->name();
    }
    if(defined($a->frequency())) {
      print " with frequency ", $a->frequency();
    }
    print "\n";



=head1 DESCRIPTION

This is a class representing a single allele of a variation.  In addition to
the nucleotide(s) (or absence of) that representing the allele frequency
and population information may be present.

=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::Allele;

use Bio::EnsEMBL::Storable;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);

our @ISA = ('Bio::EnsEMBL::Storable');


=head2 new

  Arg [-dbID]: int - unique internal identifier for the Allele
  Arg [-ADAPTOR]: Bio::EnsEMBL::Variation::DBSQL::AlleleAdaptor
  Arg [-ALLELE]: string - the nucleotide string representing the allele
  Arg [-FREQUENCY]: float - the frequency of the allele
  Arg [-POPULATION]: Bio::EnsEMBL::Variation::Population - the population
                     in which the allele was recorded
  Example    :     $allele = Bio::EnsEMBL::Variation::Allele->new
                      (-allele => 'A',
                       -frequency => 0.85,
                       -population => $pop);

  Description: Constructor.  Instantiates a new Allele object.
  Returntype : Bio::EnsEMBL::Variation::Allele
  Exceptions : none
  Caller     : general

=cut


sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my ($dbID, $adaptor, $allele, $freq, $pop) =
    rearrange(['dbID', 'ADAPTOR', 'ALLELE', 'FREQUENCY', 'POPULATION'], @_);

  return bless {'dbID'    => $dbID,
                'adaptor' => $adaptor,
                'allele'  => $allele,
                'frequency' => $freq,
                'population' => $pop}, $class;
}



=head2 allele

  Arg [1]    : string $newval (optional) 
               The new value to set the allele attribute to
  Example    : print $a->allele();
               $a1->allele('A');
               $a2->allele('-');
  Description: Getter/Setter for the allele attribute.  The allele is a string
               of nucleotide sequence, or a '-' representing the absence of
               sequence (deletion).
  Returntype : string
  Exceptions : none
  Caller     : general

=cut

sub allele{
  my $self = shift;
  return $self->{'allele'} = shift if(@_);
  return $self->{'allele'};
}




=head2 frequency

  Arg [1]    : float $newval (optional) 
               The new value to set the frequency attribute to
  Example    : $frequency = $a->frequency();
  Description: Getter/Setter for the frequency attribute. The frequency is
               the frequency of the occurance of the allele. If the population
               attribute it is the frequency of the allele within that
               population.
  Returntype : float
  Exceptions : none
  Caller     : general

=cut

sub frequency{
  my $self = shift;
  return $self->{'frequency'} = shift if(@_);
  return $self->{'frequency'};
}



=head2 population

  Arg [1]    : Bio::EnsEMBL::Variation::Population $newval (optional)
               The new value to set the population attribute to
  Example    : $population = $a->population();
  Description: Getter/Setter for the population attribute
  Returntype : Bio::EnsEMBL::Variation::Population
  Exceptions : throw on incorrect argument
  Caller     : general

=cut

sub population{
  my $self = shift;

  if(@_) {
    if(!ref($_[0]) || !$_[0]->isa('Bio::EnsEMBL::Variation::Population')) {
      throw('Bio::EnsEMBL::Variation::Population argument expected.');
    }
    $self->{'population'} = shift;
  }

  return $self->{'population'};
}


1;
