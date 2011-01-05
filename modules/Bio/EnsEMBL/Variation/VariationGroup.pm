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

# Ensembl module for Bio::EnsEMBL::Variation::VariationGroup
#
# Copyright (c) 2004 Ensembl
#

=head1 NAME

Bio::EnsEMBL::Variation::VariationGroup - Ensembl representation of a grouping of
variations (aka haplotype set).

=head1 SYNOPSIS

  use Bio::EnsEMBL::Variation::VariationGroup;

  ...



=head1 DESCRIPTION

This is a class representing a grouping of variations that have tight linkage.
This is commonly known as a Haplotype Set.  It can be viewed as a vertical
grouping of AlleleGroups.

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::Variation::VariationGroup;

use Bio::EnsEMBL::Storable;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw deprecate warning);

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Storable);



=head2 new

  Arg [dbID] :
    int - unique internal identifier for this allele group

  Arg [ADAPTOR] :
    Bio::EnsEMBL::Variation::DBSQL::VariationGroupAdaptor

  Arg [NAME] :
    string - the name of this variation group

  Arg [SOURCE] :
    string - the name of the database this variation group is from

  Arg [TYPE] :
    string - the type of this group.  Must be either 'haplotype' or 'tag'.

  Arg [VARIATIONS] :
    reference to list of Bio::EnsEMBL::Variation::Variation objects - the
    variations which make up this variation group.

  Example    :
    $ag = Bio::EnsEMBL::Variation::VariationGroup->new
      (-name   => 'PERLEGEN:B000007',
       -source => 'dbSNP',
       -type => 'haplotype',
       -variations => [$var1, $var2, $var3, $var4]);
  Description: Constructor.  Instantiates a new VariationGroup
  Returntype : Bio::EnsEMBL::Variation::VariationGroup
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub new {
  my $class = shift;

  my ($dbID, $adaptor, $name, $type, $src, $vars) =
    rearrange([qw(DBID ADAPTOR NAME TYPE SOURCE VARIATIONS)], @_);

  if(defined($vars) && !ref($vars) eq 'ARRAY') {
    throw("Reference to list of Bio::EnsEMBL::Variation::Variation ".
          "objects expected");
  }

  $vars ||= [];

  return bless {'dbID' => $dbID,
                'adaptor' => $adaptor,
                'name' => $name,
                'type' => $type,
                'source' => $src,
                'variations' => $vars}, $class;
}



=head2 name

  Arg [1]    : string $name
  Example    : print $allele_group->name();
  Description: Getter/Setter for the name of this VariationGroup
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



=head2 source

  Arg [1]    : (optional) string $source
  Example    : print $ag->source(), "\n";
  Description: Getter/Setter for the name of the source database of this
               VariationGroup.
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



=head2 type

  Arg [1]    : string $newval (optional)
               The new value to set the type attribute to
  Example    : $type = $vg->type()
  Description: Getter/Setter for the type attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub type{
  my $self = shift;
  return $self->{'type'} = shift if(@_);
  return $self->{'type'};
}


=head2 get_all_Variations

  Arg [1]    : none
  Example    :

  Description: Retrieves all variations that make up this VariationGroup.
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
  Example    : $vg->add_Variation($var2);
  Description: Adds a Variation to this VariationGroup
  Returntype : none
  Exceptions : throw on incorrect arguments
  Caller     : general
  Status     : At Risk

=cut

sub add_Variation {
  my $self = shift;
  my $var = shift;

  if(!ref($var) || !$var->isa('Bio::EnsEMBL::Variation::Variation')) {
    throw('Bio::EnsEMBL::Variation::Variation argument expected');
  }

  push @{$self->{'variations'}}, $var;

  return;
}


1;
