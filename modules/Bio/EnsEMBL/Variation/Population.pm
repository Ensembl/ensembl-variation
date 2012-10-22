=head1 LICENSE

 Copyright (c) 1999-2012 The European Bioinformatics Institute and
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

# Ensembl module for Bio::EnsEMBL::Variation::Population
#
# Copyright (c) 2004 Ensembl
#


=head1 NAME

Bio::EnsEMBL::Variation::Population - A population represents a phenotypic 
group, ethnic group, set of individuals used in an assay, etc.

=head1 SYNOPSIS

    # Population
    $pop = Bio::EnsEMBL::Variation::Population->new
       (-name => 'WEST AFRICA',
        -description => 'Sub-Saharan Nations bordering Atlantic north' .
                        ' of Congo River, and Central/Southern Atlantic' .
                        ' Island Nations.');

    ...

    # print out all sub populations of a population
    # same could work for super populations

    print_sub_pops($pop);

    sub print_sub_pops {
      my $pop = shift;
      my $level = shift || 0;

      my $sub_pops = $pop->get_all_sub_Populations();

      foreach my $sp (@$sub_pops) {
        print ' ' x $level++,
              'name: ', $sp->name(),
              'desc: ', $sp->description(),
              'size: ', $sp->size(),"\n";
        print_sub_pops($sp, $level);
      }
    }



=head1 DESCRIPTION

This is a class representing a population.  A population may consist of any
grouping of individuals, including phenotypic groups (e.g. people with
diabetes), ethnic groups (e.g. caucasians), individuals used in an assay
(e.g. subjects in experiment X), etc.

Populations may be arranged into an arbitrary hierarchy of sub and super
populations.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::Population;

use Bio::EnsEMBL::Variation::Sample;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw);

our @ISA = ('Bio::EnsEMBL::Variation::Sample');


=head2 new

  Arg [-dbID]: int - unique internal identifier of the sample
  Arg [-ADAPTOR]: Bio::EnsEMBL::PopulationAdaptor
  Arg [-NAME]: string - name of the population
  Arg [-DESCRIPTION]: string - description of the population
  Arg [-SIZE]: int - the size of the population
  Arg [-SUB_POPULATIONS]: listref of Bio::EnsEMBL::Population objects 
  Example    : $pop = Bio::EnsEMBL::Variation::Population->new
       (-name => 'WEST AFRICA',
        -description => 'Sub-Saharan Nations bordering Atlantic north' .
                        ' of Congo River, and Central/Southern Atlantic' .
                        ' Island Nations.'
        -sub_populations => \@sub_pops);
  Description: Constructor. Instantiates a new Population object
  Returntype : Bio::EnsEMBL::Variation::Population
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub new {
  my $caller = shift;

  my $class = ref($caller) || $caller;

  my ($dbID, $adaptor, $name, $desc, $size, $is_strain, $sub_pops) =
    rearrange(['DBID','ADAPTOR','NAME', 'DESCRIPTION', 'SIZE',
               'SUB_POPULATIONS'], @_);

  return bless {'dbID'        => $dbID,
                'adaptor'     => $adaptor,
                'name'        => $name,
                'description' => $desc,
                'size'        => $size,
                'sub_populations' => $sub_pops}, $class;
}



=head2 get_all_sub_Populations

  Arg [1]    : none
  Example    : foreach my $sub_pop (@{$pop->get_all_sub_Populations}) {
                 my $sub_sub_pops = $sub_pop->get_all_sub_Populations();
               }
  Description: Retrieves all populations which are conceptually a sub set
               of this population.
  Returntype : reference to list of Bio::EnsEMBL::Variation::Population objects
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_sub_Populations {
  my $self = shift;

  if(!defined($self->{'sub_populations'}) && $self->{'adaptor'}) {
    # lazy-load from database
    $self->{'sub_populations'} =
      $self->{'adaptor'}->fetch_all_by_super_Population($self);
  }
  return $self->{'sub_populations'} || [];
}



=head2 get_all_super_Populations

  Arg [1]    : none
  Example    : foreach my $sup_pop (@{$pop->get_all_super_Populations}) {
                 my $sup_sup_pops = $sup_pop->get_all_super_Populations();
               }
  Description: Retrieves all populations which this population is a part of
               from the database.
               Super populations may not be directly added in order to avoid
               circular references and memory leaks.  You must add
               sub_Populations instead and store this in the database.
  Returntype : reference to list of Bio::EnsEMBL::Variation::Population objects
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_super_Populations {
  my $self = shift;

  return [] if(!$self->{'adaptor'});

  # load from database - do not cache to avoid circular references (mem leak)!
  return $self->{'adaptor'}->fetch_all_by_sub_Population($self);
}



=head2 add_sub_Population

  Arg [1]    : Bio::EnsEMBL::Variation::Population $pop
  Example    : $pop->add_sub_Population($sub_pop);
               $sub_pop->add_super_Population($pop);
  Description: Adds a sub population to this population.
  Returntype : none
  Exceptions : throw on incorrect argument
  Caller     : general
  Status     : At Risk

=cut

sub add_sub_Population {
  my $self = shift;
  my $pop = shift;

  if(!ref($pop) || !$pop->isa('Bio::EnsEMBL::Variation::Population')) {
    throw('Bio::EnsEMBL::Variation::Population argument expected.');
  }

  if($pop == $self) {
    throw("Cannot add self as sub population.");
  }

  $self->{'sub_populations'} ||= [];
  push @{$self->{'sub_populations'}}, $pop;

  return $pop;
}

=head2 get_all_synonyms

  Arg [1]    : (optional) string $source - the source of the synonyms to
               return.
  Example    : @dbsnp_syns = @{$p->get_all_synonyms('dbSNP')};
               @all_syns = @{$p->get_all_synonyms()};
  Description: Retrieves synonyms for this Population. If a source argument
               is provided all synonyms from that source are returned,
               otherwise all synonyms are returned.
  Returntype : reference to list of strings
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub get_all_synonyms {
  my $self = shift;
  my $source = shift;

  return [] if(!$self->adaptor()); #if there is no adaptor, return empty string

  return $self->adaptor()->fetch_synonyms($self->dbID(),$source);

}

=head2 get_all_Individuals

  Arg [1]    : none
  Example    : @individuals = @{$p->get_all_individuals()};
  Description: Retrieves all Individuals belonging to this Population.
  Returntype : reference to list of Bio::EnsEMBL::Variation::Individual objects
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub get_all_Individuals {
  my $self = shift;

  my $ia = $self->adaptor->db->get_IndividualAdaptor;
  
  return (defined $ia ? $ia->fetch_all_by_Population($self) : []);
}

1;
