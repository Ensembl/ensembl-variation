=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut


=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

 Questions may also be sent to the Ensembl help desk at
 <http://www.ensembl.org/Help/Contact>.

=cut

# Ensembl module for Bio::EnsEMBL::Variation::Allele
#
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

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::Allele;

use Bio::EnsEMBL::Storable;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Scalar qw(assert_ref check_ref);
use Scalar::Util qw(weaken);
use Bio::EnsEMBL::Variation::Failable;

our @ISA = ('Bio::EnsEMBL::Storable', 'Bio::EnsEMBL::Variation::Failable');


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
  Status     : Stable

=cut


sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my ($dbID, $adaptor, $allele, $freq, $count, $pop, $ss_id, $variation_id, $population_id) =
    rearrange(['dbID', 'ADAPTOR', 'ALLELE', 'FREQUENCY', 'COUNT', 'POPULATION', 'SUBSNP', 'VARIATION_ID', 'POPULATION_ID'], @_);
  
  # set subsnp_id to undefined if it's 0 in the DB
  #$ss_id = undef if (defined $ss_id && $ss_id == 0);
  
  # add ss to the subsnp_id
  $ss_id = 'ss'.$ss_id if defined $ss_id && $ss_id !~ /^ss/;

  # Check that we at least get a BaseAdaptor
  assert_ref($adaptor,'Bio::EnsEMBL::Variation::DBSQL::BaseAdaptor');
  # If the adaptor is not an AlleleAdaptor, try to get it via the passed adaptor
  unless (check_ref($adaptor,'Bio::EnsEMBL::Variation::DBSQL::AlleleAdaptor')) {
      $adaptor = $adaptor->db->get_AlleleAdaptor();
      # Verify that we could get the AlleleAdaptor
        assert_ref($adaptor,'Bio::EnsEMBL::Variation::DBSQL::AlleleAdaptor');
  }
  
  my $self = bless {}, $class;
  
  $self->dbID($dbID);
  $self->adaptor($adaptor);
  $self->allele($allele);
  $self->frequency($freq);
  $self->count($count);
  $self->subsnp($ss_id);
  $self->{'_variation_id'} = $variation_id;
  $self->{'_population_id'} = $population_id;
  $self->population($pop) if (defined($pop));
  
  return $self;
}

sub new_fast {
  my $class = shift;
  my $hashref = shift;
  return bless $hashref, $class;
}


# An internal method for getting a unique hash key identifier, used by the Variation module 
sub _hash_key {
    my $self = shift;
    
    # By default, return the dbID
    my $dbID = $self->dbID();
    return $dbID if (defined($dbID));
     
    # If no dbID is specified, e.g. if we are creating a 'custom' object, return a fake dbID. This is necessary since e.g. Variation stores
    # its alleles in a hash with dbID as key. To create fake dbIDs, use the string representing the memory address.
    ($dbID) = sprintf('%s',$self) =~ m/\(([0-9a-fx]+)\)/i;
    return $dbID;
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
  Status     : Stable

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
  Status     : Stable

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
  Status     : Stable

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
  Status     : Stable

=cut

sub population{
    my $self = shift;

    if(@_) {
        assert_ref($_[0],'Bio::EnsEMBL::Variation::Population');
        $self->{'population'} = shift;
        $self->{'_population_id'} = $self->{'population'}->dbID();
    }

    # Population can be lazy-loaded, so get it from the database if we have a population_id but no cached object
    if (!defined($self->{'population'}) && defined($self->{'_population_id'})) {
        
        # Check that an adaptor is attached
        assert_ref($self->adaptor(),'Bio::EnsEMBL::Variation::DBSQL::AlleleAdaptor');
        
        # Get a population object
        my $population = $self->adaptor->db->get_PopulationAdaptor()->fetch_by_dbID($self->{'_population_id'});
        
        # Set the population
				$self->{'population'} = $population;
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
  Status     : Stable

=cut

sub subsnp{
  my $self = shift;
  if(@_) {
    $self->{'subsnp'} = shift;
  }
  
  my $ssid = $self->{'subsnp'};
  if(defined($ssid)) {
	$ssid = 'ss'.$ssid unless $ssid =~ /^ss/;
  }
  
  return $ssid;
}


=head2 variation

  Arg [1]    : Bio::EnsEMBL::Variation::Variation $newval (optional) 
               The new value to set the variation attribute to
  Example    : print $a->variation->name();
  Description: Getter/Setter for the variation attribute.
  Returntype : Bio::EnsEMBL::Variation::Variation
  Exceptions : throw on incorrect argument
  Caller     : general

=cut

sub variation {
    my $self = shift;
    my $variation = shift;
  
    # Set the dbID of the variation object on this allele
    if(defined($variation)) {
        assert_ref($variation,'Bio::EnsEMBL::Variation::Variation');
        $self->{'_variation_id'} = $variation->dbID();
        return $variation;
    }

    # Load the variation from the database if we have a variation_id
    if (defined($self->{'_variation_id'})) {
        
        # Check that an adaptor is attached
        assert_ref($self->adaptor(),'Bio::EnsEMBL::Variation::DBSQL::BaseAdaptor');
        
        # Get a variation object
        $variation = $self->adaptor->db->get_VariationAdaptor()->fetch_by_dbID($self->{'_variation_id'});
		
		$self->{variation} = $variation;
    }
    
    # Return the variation object
    return $self->{variation};
}


=head2 subsnp_handle

  Arg [1]    : string $newval (optional) 
               The new value to set the subsnp_handle attribute to
  Example    : print $a->subsnp_handle();
  Description: Getter/Setter for the subsnp_handle attribute.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub subsnp_handle{
    my $self = shift;
    my $handle = shift;
      
    # if changing handle
    if(defined($handle)) {
        $self->{'subsnp_handle'} = $handle;
    }
    elsif (!exists($self->{'subsnp_handle'})) {

        # Check that this allele has an adaptor attached
        assert_ref($self->adaptor(),'Bio::EnsEMBL::Variation::DBSQL::AlleleAdaptor');
        
        $self->{'subsnp_handle'} = $self->adaptor->get_subsnp_handle($self);
    }
    
    return $self->{'subsnp_handle'};
}

=head2 frequency_subsnp_handle

  Arg[1]     : Bio::EnsEMBL::Variation::Population
                       The population object to get/set the frequency submitter subsnp handle for      
  Arg[2]     : string $newval (optional) 
               The new value to set the frequency_subsnp_handle attribute to
  Example    : print $a->frequency_subsnp_handle($population);
  Description: Getter/Setter for the frequency_subsnp_handle attribute.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub frequency_subsnp_handle{
    my $self = shift;
    my $population = shift;
    my $handle = shift;

    assert_ref($population,'Bio::EnsEMBL::Variation::Population') if (defined($population));
      
    # if changing handle
    if(defined($handle)) {
        $self->{'frequency_subsnp_handle'} = $handle;
    }
    elsif (!exists($self->{'frequency_subsnp_handle'})) {

        # Check that this allele has an adaptor attached
        assert_ref($self->adaptor(),'Bio::EnsEMBL::Variation::DBSQL::AlleleAdaptor');
        
        $self->{'frequency_subsnp_handle'} = $self->adaptor->get_subsnp_handle($self,$population );
    }
    
    return $self->{'frequency_subsnp_handle'};
}

sub _weaken {
    my $self = shift;
    
    # If the variation is not defined, do nothing
    return unless (defined($self->variation()));
    
    # Weaken the link to the variation
    weaken($self->{'variation'});
}

1;
