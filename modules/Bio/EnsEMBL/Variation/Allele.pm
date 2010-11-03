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
  Status     : At Risk

=cut


sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my ($dbID, $adaptor, $allele, $freq, $count, $pop, $ss_id) =
    rearrange(['dbID', 'ADAPTOR', 'ALLELE', 'FREQUENCY', 'COUNT', 'POPULATION', 'SUBSNP'], @_);
  
  # set subsnp_id to undefined if it's 0 in the DB
  $ss_id = undef if (defined $ss_id && $ss_id == 0);
  
  # add ss to the subsnp_id
  $ss_id = 'ss'.$ss_id if defined $ss_id && $ss_id !~ /^ss/;

  return bless {'dbID'    => $dbID,
                'adaptor' => $adaptor,
                'allele'  => $allele,
                'frequency' => $freq,
                'count'   => $count,
                'population' => $pop,
                'subsnp'  => $ss_id}, $class;
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
  Status     : At Risk

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
  Example    : $frequency = $allele->count()
  Description: Getter/Setter for the observed count of this allele
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



=head2 population

  Arg [1]    : Bio::EnsEMBL::Variation::Population $newval (optional)
               The new value to set the population attribute to
  Example    : $population = $a->population();
  Description: Getter/Setter for the population attribute
  Returntype : Bio::EnsEMBL::Variation::Population
  Exceptions : throw on incorrect argument
  Caller     : general
  Status     : At Risk

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
      
      my $row = $sth->fetchrow_arrayref();
      
      my $handle;
      
      $handle = $row->[0] if defined($row);
      
      return $self->{'subsnp_handle'} = $handle;
    }
  }
  
  return $self->{'subsnp_handle'};
}


1;
