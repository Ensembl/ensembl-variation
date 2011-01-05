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

# Ensembl module for Bio::EnsEMBL::Variation::AlleleGoup
#
# Copyright (c) 2004 Ensembl
#

=head1 NAME

Bio::EnsEMBL::Variation::AlleleGroup - Ensembl representation of a grouping of
alleles (aka haplotype).

=head1 SYNOPSIS

  use Bio::EnsEMBL::Variation::AlleleGroup;

  ...

  # instantiate a new AlleleGroup
  $ag = Bio::EnsEMBL::Variation::AlleleGroup->new
    (-name   => 'ABDR-20',
     -variation_group => $var_group,
     -population => $pop,
     -source => 'dbSNP',
     -frequency => 0.32);


    # add some variations and the alleles that make up this group

    $ag->add_Variation($var1, 'T');
    $ag->add_Variation($var2, 'C');

   # print some info about an allele group

    print $ag->name(), "\n";
    print $ag->variation_group()->name(), "\n";
    print $ag->population()->name(), "\n" if($ag->population());
    print $ag->source(), "\n";
    print $ag->frequency(), "\n" if($ag->frequency());

    # print out the names of the variations and the alleles that comprise
    # this allele group

    $vars = $ag->get_all_Variations();
    $alleles = $ag->get_all_alleles();

    for($i=0; $i < @$vars; $i++) {
      print $vars->[$i]->name(), ':', $alleles->[$i], "\n";
    }


=head1 DESCRIPTION

This is a class representing a grouping of alleles that have tight linkage and
are usually present together.  This is commonly known as a Haplotype or
Haplotype Block.

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::Variation::AlleleGroup;

use Bio::EnsEMBL::Storable;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw deprecate warning);

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Storable);



=head2 new

  Arg [dbID] :
    int - unique internal identifier for this allele group

  Arg [ADAPTOR] :
    Bio::EnsEMBL::Variation::DBSQL::AlleleGroupAdaptor

  Arg [NAME] :
    string - the name of this allele group

  Arg [POPULATION] :
    Bio::EnsEMBL::Variation::Population - the population this AlleleGroup is
    found in.

  Arg [SOURCE] :
    string - the name of the database this allele group is from

  Arg [VARIATION_GROUP] :
    Bio::EnsEMBL::Variation::VariationGroup - the horizontal grouping of
    variations that this allele group belongs to.

  Arg [FREQUENCY] :
    The frequency of occurance of this AlleleGroup

  Example    :
    $ag = Bio::EnsEMBL::Variation::AlleleGroup->new
      (-name   => 'ABDR-20',
       -variation_group => $var_group,
       -population => $pop,
       -source => 'dbSNP',
       -frequency => 0.32);
  Description: Constructor.  Instantiates a new AlleleGroup
  Returntype : Bio::EnsEMBL::Variation::AlleleGroup
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub new {
  my $class = shift;

  my ($dbID, $adaptor, $name, $pop, $src, $var_grp, $freq) =
    rearrange([qw(DBID ADAPTOR NAME POPULATION 
                  SOURCE VARIATION_GROUP FREQUENCY)], @_);

  if(defined($pop) && 
     (!ref($pop) || !$pop->isa('Bio::EnsEMBL::Variation::Population'))) {
    throw('Bio::EnsEMBL::Variation::Population argument expected');
  }

  if(defined($var_grp) && (!ref($var_grp) || 
     !$var_grp->isa('Bio::EnsEMBL::Variation::VariationGroup'))) {
    throw('Bio::EnsEMBL::Variation::VariationGroup argument expected');
  }

  return bless {'dbID' => $dbID,
                'adaptor' => $adaptor,
                'name' => $name,
                'population' => $pop,
                'source' => $src,
                'variation_group' => $var_grp,
                'frequency' => $freq}, $class;
}



=head2 name

  Arg [1]    : string $name
  Example    : print $allele_group->name();
  Description: Getter/Setter for the name of this AlleleGroup
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub name {
  my $self = shift;
  $self->{'name'} = shift if(@_);
  return $self->{'name'};
}



=head2 population

  Arg [1]    : (optional) Bio::EnsEMBL::Variation::Population $pop
  Example    : print $ag->population()->name() if($ag->population());
  Description: Getter/Setter for the Population associated with this
               allele group.  Returns undef if there is no population
               associated with this allele group
  Returntype : Bio::EnsEMBL::Variation::Population
  Exceptions : tr
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
    $self->{'population'} = $pop;
  }

  return $self->{'population'};
}



=head2 variation_group

  Arg [1]    : (optional) Bio::EnsEMBL::Variation::VariationGroup $vg
  Example    : print $ag->variation_group()->name();
  Description: Getter/Setter for VariationGroup associated with this allele
               group.
  Returntype : Bio::EnsEMBL::Variation::VariationGroup
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub variation_group {
  my $self = shift;
  if(@_) {
    my $vg = shift;
    if(defined($vg) &&
       (!ref($vg) || !$vg->isa('Bio::EnsEMBL::Variation::VariationGroup'))) {
      throw('Bio::EnsEMBL::Variation::VariationGroup argument expected');
    }
    $self->{'variation_group'} = $vg
  }
  return $self->{'variation_group'};
}



=head2 source

  Arg [1]    : (optional) string $source
  Example    : print $ag->source(), "\n";
  Description: Getter/Setter for the name of the source database of this
               AlleleGroup.
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub source {
  my $self = shift;
  $self->{'source'} = shift if(@_);
  return $self->{'source'};
}


=head2 frequency

  Arg [1]    : (optional) float $freq
  Example    : print $ag->frequency(), "\n";
  Description: Getter/Setter for the frequence of this AlleleGroup
  Returntype : float
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub frequency {
  my $self = shift;
  $self->{'frequency'} = shift if(@_);
  return $self->{'frequency'};
}




=head2 get_all_alleles

  Arg [1]    : none
  Example    :
    $vars = $ag->get_all_Variations();
    $alleles = $ag->get_all_alleles();

    for($i=0; $i < @$vars; $i++) {
      print $vars->[$i]->name(), ':', $alleles->[$i], "\n";
    }
  Description: Retrieves all alleles that are part of this allele group.
               The alleles are returned as simple strings.  The ordering of
               the alleles returned by this method are consistant with the
               ordering of the variations returned by the get_all_Variations
               method.  The third allele corresponds to the third variation,
               etc.
  Returntype : reference to a list of strings
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub get_all_alleles {
  my $self = shift;
  return $self->{'alleles'};
}



=head2 get_all_Variations

  Arg [1]    : none
  Example    :
    $vars = $ag->get_all_Variations();
    $alleles = $ag->get_all_alleles();

    for($i=0; $i < @$vars; $i++) {
      print $vars->[$i]->name(), ':', $alleles->[$i], "\n";
    }
  Description: Retrieves all variations that are associated with the alleles
               that comprise this allele group.
               The ordering of the Variations returned by this method are
               consistant with the ordering of the alleles returned by the
               get_all_alleles method.  The third allele corresponds to the
               third variation, etc.
  Returntype : reference to a list of strings
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub get_all_Variations {
  my $self = shift;
  return $self->{'variations'};
}



=head2 add_Variation

  Arg [1]    : Bio::EnsEMBL::Variation::Variation $var
  Arg [2]    : string $allele
  Example    : $ag->add_Variation($var2, 'C');
  Description: Adds an allele and associated variation to this allele group
  Returntype : none
  Exceptions : throw on incorrect arguments
  Caller     : AlleleGroupAdaptor
  Status     : At Risk

=cut

sub add_Variation {
  my $self = shift;
  my $var = shift;
  my $allele = shift;

  if(!ref($var) || !$var->isa('Bio::EnsEMBL::Variation::Variation')) {
    throw('Bio::EnsEMBL::Variation::Variation argument expected');
  }

  if(!$allele) {
    throw('allele argument expeceted');
  }

  push @{$self->{'variations'}}, $var;
  push @{$self->{'alleles'}}, $allele;

  return;
}


1;
