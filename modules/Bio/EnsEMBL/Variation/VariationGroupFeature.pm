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

# Ensembl module for Bio::EnsEMBL::Variation::VariationGroupFeature
#
# Copyright (c) 2004 Ensembl
#


=head1 NAME

Bio::EnsEMBL::Variation::VariationGroupFeature - A genomic position for a variation group (aka haplotype block).

=head1 SYNOPSIS

    # VariationGroup feature representing the location of a haplotype block
    $vgf = Bio::EnsEMBL::Variation::VariationGroupFeature->new
       (-start   => 100,
        -end     => 100,
        -strand  => 1,
        -slice   => $slice,
        -variation_group => $vg);

    ...

    # a variation group feature is like any other ensembl feature, can be
    # transformed etc.
    $vgf = $vgf->transform('supercontig');

    print $vgf->start(), "-", $vgf->end(), '(', $vgf->strand(), ')', "\n";

    print $vgf->variation_group()->name(), "\n";

=head1 DESCRIPTION

This is a class representing the genomic position of a VariationGroup.  A
VariationGroup is a collection of tightly linked Variations sometimes
known as a Haplotype block.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::VariationGroupFeature;

use Bio::EnsEMBL::Feature;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument  qw(rearrange);


our @ISA = ('Bio::EnsEMBL::Feature');


=head2 new

  Arg [-dbID] :
    see superclass constructor

  Arg [-ADAPTOR] :
    see superclass constructor

  Arg [-START] :
    see superclass constructor
  Arg [-END] :
    see superclass constructor

  Arg [-STRAND] :
    see superclass constructor

  Arg [-SLICE] :
    see superclass constructor

  Arg [-VARIATION_GROUP] :
    Bio::EnsEMBL::Variation::VariationGroup - the variation group that this
    feature represents the genomic position of.

  Arg [-VARIATION_GROUP_NAME] :
    string - the name of this variation group.  This allows the name to be
    rapidly obtained even when the VariationGroup has not been retrieved from
    the database.

  Arg [-VARIATION_GROUP_ID] :
    int - the internal identifier of the variation group.  This can be set
    instead of the VariationGroup to allow the variation group to be
    lazy-loaded later.

  Example    :
    $vgf = Bio::EnsEMBL::Variation::VariationGroupFeature->new
       (-start   => 100,
        -end     => 100,
        -strand  => 1,
        -slice   => $slice,
        -variation_group_name => $vg->name(),
        -variation_group => $vg);


  Description: Constructor. Instantiates a new VariationFeature object.
  Returntype : Bio::EnsEMBL::Variation::Variation
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);

  my ($vg, $vg_name, $vg_id) =
   rearrange([qw(VARIATION_GROUP VARIATION_GROUP_NAME VARIATION_GROUP_ID)],@_);

  $self->{'variation_group'} = $vg;
  $self->{'variation_group_name'}  = $vg_name;
  $self->{'_variation_group_id'} = $vg_id;

  return $self;
}


=head2 variation_group_name

  Arg [1]    : string $newval (optional) 
               The new value to set the variation_group_name attribute to
  Example    : $variation_group_name = $obj->variation_group_name()
  Description: Getter/Setter for the variation_group_name attribute
  Returntype : string
  Exceptions : none
  Caller     : general
  Status     : At Risk

=cut

sub variation_group_name{
  my $self = shift;
  return $self->{'variation_group_name'} = shift if(@_);
  return $self->{'variation_group_name'};
}


=head2 variation_group

  Arg [1]    : Bio::EnsEMBL::Variation::VaritionGroup $vg (optional)
  Example    : print $vg->variation_group->name();
  Description: Getter/Setter for the VariationGroup associated with this
               feature.  If not set and this feature has an associated adaptor,
               an attempt will be made to lazy-load the variation from the
               database.
  Returntype : Bio::EnsEMBL::Variation::VariationGroup
  Exceptions : throw on bad argument
  Caller     : general
  Status     : At Risk

=cut

sub variation_group {
  my $self = shift;


  if(@_) {
    if(!ref($_[0]) || !$_[0]->isa('Bio::EnsEMBL::Variation::VariationGroup')) {
      throw("Bio::EnsEMBL::Variation::VariationGroup argument expected");
    }
    $self->{'variation_group'} = shift;
  }
  elsif(!defined($self->{'variation_group'}) && $self->{'adaptor'} &&
        defined($self->{'_variation_group_id'})) {
    # lazy-load from database on demand
    my $vga = $self->{'adaptor'}->db()->get_VariationGroupAdaptor();
    $self->{'variation_group'} = 
      $vga->fetch_by_dbID($self->{'_variation_group_id'});
  }

  return $self->{'variation_group'};
}


=head2 display_id

  Arg [1]    : none
  Example    : print $vgf->display_id(), "\n";
  Description: Returns the 'display' identifier for this feature. For
               VariationGroupFeatures this is simply the name of the variation
               group it is associated with.
  Returntype : string
  Exceptions : none
  Caller     : webcode
  Status     : At Risk

=cut

sub display_id {
  my $self = shift;
  return $self->{'variation_group_name'} || '';
}



1;
