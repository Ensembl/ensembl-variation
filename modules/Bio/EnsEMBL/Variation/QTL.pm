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

# Ensembl module for Bio::EnsEMBL::Variation::Phenotype
#
# Copyright (c) 2011 Ensembl
#


=head1 NAME

Bio::EnsEMBL::Variation::QTL - Ensembl representation of a quantitative trait
locus (QTL).

=head1 SYNOPSIS

    my $qtl = Bio::EnsEMBL::Variation::QTL->new(-stable_id => 'QTL00123');

=head1 DESCRIPTION

This is a class representing a QTL.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::QTL;

use Bio::EnsEMBL::Storable;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument  qw(rearrange);

our @ISA = ('Bio::EnsEMBL::Storable');


=head2 new
  Arg [-stable_id] : stable identifier or name for QTL
  Example : my $qtl = Bio::EnsEMBL::Variation::QTL->new(-stable_id => 'QTL00123');

  Description: Constructor. Instantiates a new QTL object.
  Returntype : Bio::EnsEMBL::Variation::QTL
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub new {
    my $caller = shift;
    my $class  = ref($caller) || $caller;
    my $self = $class->SUPER::new(@_);
    my ($stable_id) = rearrange([qw(STABLE_ID)], @_);
    $self = {
        'stable_id' => $stable_id,
    };
    return bless $self, $class;
}

sub new_fast {
  my $class = shift;
  my $hashref = shift;
  return bless $hashref, $class;
}

sub stable_id {
    my $self = shift;
    return $self->{stable_id} = shift if @_;
    return $self->{stable_id};
}

1;