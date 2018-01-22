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

# Ensembl module for Bio::EnsEMBL::Variation::IndividualGenotype
#
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

package Bio::EnsEMBL::Variation::IndividualGenotype;

use Bio::EnsEMBL::Variation::Genotype;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw deprecate warning);

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Variation::Genotype);

=head2 new

  Arg [-adaptor] :
    Bio::EnsEMBL::Variation::DBSQL::IndividualGenotypeAdaptor
  Arg [-genotype] :
    arrayref - arrayref of alleles making up this genotype (in haplotype order)
  Arg [-allele2] :
    string - One of the two alleles defining this genotype
  Arg [-variation] :
    Bio::EnsEMBL::Variation::Variation - The variation associated with this
    genotype
  Arg [-individual] :
    Bio::EnsEMBL::Individual - The individual this genotype is for.
  Example    : $ind_genotype = Bio:EnsEMBL::Variation::IndividualGenotype->new(
				-start      => 100,
				-end        => 100,
				-strand     => 1,
				-slice      => $slice,
				-genotype   => ['A','T'],
				-variation  => $variation,
				-individual => $ind
			   );
  Description: Constructor.  Instantiates an IndividualGenotype object.
  Returntype : Bio::EnsEMBL::Variation::IndividualGenotype
  Exceptions : throw on bad argument
  Caller     : general
  Status     : Deprecated. Use Bio::EnsEMBL::Variation::SampleGenotype::new.

=cut

sub new {
 my $caller = shift;
  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);

  my ($adaptor, $genotype, $var, $varid, $ssid, $ind, $phased) =
    rearrange([qw(adaptor genotype variation _variation_id subsnp individual phased)],@_);

  if(defined($var) && (!ref($var) || !$var->isa('Bio::EnsEMBL::Variation::Variation'))) {
    throw("Bio::EnsEMBL::Variation::Variation argument expected");
  }

  if(defined($ind) && (!ref($ind) || !$ind->isa('Bio::EnsEMBL::Variation::Individual'))) {
    throw("Bio::EnsEMBL::Variation::Individual argument expected");
  }

  $self->{'adaptor'}       = $adaptor;
  $self->{'genotype'}      = $genotype;
  $self->{'individual'}    = $ind;
  $self->{'variation'}     = $var;
  $self->{'subsnp'}        = $ssid;
  $self->{'phased'}        = $phased;
  $self->{'_variation_id'} = $varid unless defined $var;

  deprecate("Please use Bio::EnsEMBL::Variation::SampleGenotype::new.\n");

  return $self;
}


=head2 individual

  Arg [1]    : (optional) Bio::EnsEMBL::Variation::Individual $ind
  Example    : $ind = $ind_genotype->individual();
  Description: Getter/Setter for the individual associated with this genotype
  Returntype : Bio::EnsEMBL::Variation::Individual
  Exceptions : throw on bad argument
  Caller     : general
  Status     : Deprecated. Use Bio::EnsEMBL::Variation::SampleGenotype::sample.

=cut

sub individual {

  my $self = shift;
  if(@_) {
    my $ind = shift;
    if(defined($ind) && (!ref($ind) || !$ind->isa('Bio::EnsEMBL::Variation::Individual'))) {
      throw('Bio::EnsEMBL::Variation::Individual argument expected');
    }
    return $self->{'individual'} = $ind;
  }

  if(!defined($self->{individual}) && defined($self->{sample_id})) {
    my $ia = $self->adaptor->db->get_IndividualAdaptor;

    if(defined($ia)) {
      my $i = $ia->fetch_by_dbID($self->{sample_id});

      if(defined($i)) {
        $self->{individual} = $i;
        delete $self->{sample_id};
      }
    }
  }
  deprecate("Please use Bio::EnsEMBL::Variation::SampleGenotype::sample.\n");

  return $self->{'individual'};
}

1;
