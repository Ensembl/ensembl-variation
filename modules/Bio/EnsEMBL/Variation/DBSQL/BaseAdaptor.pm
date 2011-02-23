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


# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::BaseAdaptor
#
# Copyright (c) 2010 Ensembl
#
# You may distribute this module under the same terms as perl itself
#
#

=head1 NAME

Bio::EnsEMBL::Variation::DBSQL::BaseAdaptor

=head1 SYNOPSIS

Abstract class - should not be instantiated.  Implementation of
abstract methods must be performed by subclasses.

=head1 DESCRIPTION

This adaptor provides generic database connectivity for various Variation objects.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::BaseAdaptor;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;

our @ISA = ('Bio::EnsEMBL::DBSQL::BaseAdaptor');

=head2 fetch_all

  Description: Returns a listref of all germline variation features
  Returntype : listref of VariationFeatures
  Status     : At risk

=cut

sub fetch_all {
    my $self = shift;
    my $constraint = 's.somatic = 0';
    return $self->generic_fetch($constraint);
}

=head2 fetch_all_somatic

  Description: Returns a listref of all somatic variation features
  Returntype : listref of VariationFeatures
  Status     : At risk

=cut

sub fetch_all_somatic {
    my $self = shift;
    my $constraint = 's.somatic = 1';
    return $self->generic_fetch($constraint);
}

# returns a hash mapping SO accessions to SO and ensembl display terms
sub _SO_mappings {
    my $self = shift;
    
    unless ($self->{_SO_mappings}) {
        $self->{_SO_mappings} = $self->db->get_AttributeAdaptor->fetch_SO_mappings;
    }
    
    return $self->{_SO_mappings};
}

sub _display_term_for_SO_accession {
    my ($self, $SO_accession, $is_somatic) = @_;
   
    my $term = $self->_SO_mappings->{$SO_accession}->{display_term};
   
    if ($is_somatic) {
        $term = 'SNV' if $term eq 'SNP';
        $term = 'somatic_'.$term;
    }
    
    return $term;
}

sub _SO_term_for_SO_accession {
    my ($self, $SO_accession, $is_somatic) = @_;
    return $self->_SO_mappings->{$SO_accession}->{SO_term}
}

sub _display_term_for_SO_term {
    my ($self, $SO_term, $is_somatic) = @_;
    my $SO_accession = $self->_SO_mappings->{$SO_term}->{SO_accession};
    return $self->_display_term_for_SO_accession($SO_accession, $is_somatic);
}

1;

