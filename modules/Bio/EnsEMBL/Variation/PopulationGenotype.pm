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

# Ensembl module for Bio::EnsEMBL::Variation::PopulationGenotype
#
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
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Variation::Genotype);



=head2 new

  Arg [-dbID] :
    int - unique internal identifier
  Arg [-adaptor] :
    Bio::EnsEMBL::Variation::DBSQL::PopulationAdaptor
  Arg [-genotype] :
    arrayref of strings - The alleles defining this genotype
  Arg [-variation] :
    Bio::EnsEMBL::Variation::Variation - The variation associated with this
    genotype
  Arg [-population] :
    Bio::EnsEMBL::Population - The population this genotype is for.
  Arg [-frequency] :
    int - the frequency this genotype occurs in this population
  Example    : $pop_genotype = Bio:EnsEMBL::Variation::PopulationGenotype->new
                   (-genotype => ['A','T'],
                    -variation => $variation,
                    -population => $pop
                    -frequency  => 0.87);
  Description: Constructor.  Instantiates a PopulationGenotype object.
  Returntype : Bio::EnsEMBL::Variation::PopulationGenotype
  Exceptions : throw on bad argument
  Caller     : general
  Status     : Stable

=cut

sub new {
  my $class = shift;

  my ($dbID, $adaptor, $genotype, $var, $pop, $freq, $count, $var_id, $ss_id) =
    rearrange([qw(dbID adaptor genotype variation population frequency count _variation_id subsnp)],@_);

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

  return bless {
    'dbID'          => $dbID,
    'adaptor'       => $adaptor,
    'genotype'      => $genotype,
    'variation'     => $var,
    '_variation_id' => defined($var) ? undef : $var_id,
    'population'    => $pop,
    'frequency'     => $freq,
    'count'         => $count,
    'subsnp'        => $ss_id
  }, $class;
}




=head2 population

  Arg [1]    : (optional) Bio::EnsEMBL::Variation::Population $pop
  Example    : $pop = $pop_genotype->population();
  Description: Getter/Setter for the population associated with this genotype
  Returntype : Bio::EnsEMBL::Variation::Population
  Exceptions : throw on bad argument
  Caller     : general
  Status     : Stable

=cut


sub population {
  my $self = shift;
  if (@_) {
    my $pop = shift;
    if (defined($pop) &&
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
  Status     : Stable

=cut

sub frequency {
  my $self = shift;
  return $self->{'frequency'} = shift if (@_);
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
  Status     : Stable

=cut

sub count {
  my $self = shift;
  return $self->{'count'} = shift if (@_);
  return $self->{'count'};
}

1;
