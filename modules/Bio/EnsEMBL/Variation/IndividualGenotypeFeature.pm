=head1 LICENSE

 Copyright (c) 1999-2013 The European Bioinformatics Institute and
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

use Bio::EnsEMBL::Variation::IndividualGenotype;
use Bio::EnsEMBL::Feature;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw deprecate warning);

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Variation::IndividualGenotype Bio::EnsEMBL::Feature);



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
  Status     : Stable

=cut

sub new {
	my $caller = shift;
	my $class = ref($caller) || $caller;
	
	my $self = $class->SUPER::new(@_);
	
	my ($adaptor, $allele1, $allele2, $var, $var_id, $ind) =
	rearrange([qw(adaptor genotype variation _variation_id individual)],@_);
	
	if(defined($var) && (!ref($var) || !$var->isa('Bio::EnsEMBL::Variation::Variation'))) {
		throw("Bio::EnsEMBL::Variation::Variation argument expected");
	}
	
	if(defined($ind) && (!ref($ind) || !$ind->isa('Bio::EnsEMBL::Variation::Individual'))) {
		throw("Bio::EnsEMBL::Variation::Individual argument expected");
	}
	
	$self->{'adaptor'} = $adaptor;
	$self->{'allele1'} = $allele1;
	$self->{'allele2'} = $allele2;
	$self->{'individual'} = $ind;
	$self->{'variation'} = $var;
	$self->{'_variation_id'} = $var_id;
	
	return $self;
}

1;
