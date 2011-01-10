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


# Ensembl module for Bio::EnsEMBL::Variation::DBSQL::BaseFeatureAdaptor
#
# Copyright (c) 2010 Ensembl
#
# You may distribute this module under the same terms as perl itself
#
#

=head1 NAME

Bio::EnsEMBL::Variation::DBSQL::BaseVariationAdaptor

=head1 SYNOPSIS

Abstract class - should not be instantiated.  Implementation of
abstract methods must be performed by subclasses.

=head1 DESCRIPTION

This adaptor provides generic database connectivity for various Variation objects.

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::BaseVariationAdaptor;

use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;

our @ISA = ('Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor');

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

# returns a hash mapping SO ids to SO and ensembl display terms
sub _variation_classes {
    my $self = shift;
    
    unless ($self->{_variation_classes}) {
        
        # fetch the mapping from the database
        
        my $dbh = $self->dbc->db_handle;
        
        my $sth = $dbh->prepare(qq{
            SELECT so_id, so_term, display_term FROM variation_class
        });
        
        $sth->execute;
        
        while (my ($SO_id, $display_term, $SO_term) = $sth->fetchrow_array) {
            $self->{_variation_classes}->{$SO_id}->{SO_term}        = $SO_term;
            $self->{_variation_classes}->{$SO_id}->{display_term}   = $display_term;
            $self->{_variation_classes}->{$SO_term}->{SO_id}        = $SO_id;
        }
    }
    
    return $self->{_variation_classes};
}

sub _display_term_for_SO_id {
    my ($self, $SO_id, $is_somatic) = @_;
    
    my $term = $self->_variation_classes->{$SO_id}->{display_term};
    
    if ($is_somatic) {
        $term = 'SNV' if $term eq 'SNP';
        $term = 'somatic_'.$term;
    }
    
    return $term;
}

sub _SO_term_for_SO_id {
    my ($self, $SO_id, $is_somatic) = @_;
    return $self->_variation_classes->{$SO_id}->{SO_term}
}

sub _display_term_for_SO_term {
    my ($self, $SO_term, $is_somatic) = @_;
    my $SO_id = $self->_variation_classes->{$SO_term}->{SO_id};
    return $self->_display_term_for_SO_id($SO_id, $is_somatic);
}

1;

