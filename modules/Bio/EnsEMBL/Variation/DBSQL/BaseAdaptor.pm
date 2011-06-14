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

use Bio::EnsEMBL::Variation::Utils::Constants qw(%OVERLAP_CONSEQUENCES);

our @ISA = ('Bio::EnsEMBL::DBSQL::BaseAdaptor');

sub AttributeAdaptor {
    my $self = shift;
    
    unless ($self->{_attribute_adaptor}) {
        $self->{_attribute_adaptor} = $self->db->get_AttributeAdaptor if $self->db;
    }
    
    return $self->{_attribute_adaptor};
}

sub _consequence_type_map {
    
    # return a hash mapping between the string terms of a mysql set and
    # the equivalent numerical values

    my ($self, $table, $column) = @_;
    
    my $key = $table.'_'.$column.'_map';

    unless ($self->{$key}) {

        my $map_sth = $self->prepare(qq{SHOW COLUMNS FROM $table LIKE '$column'});
        
        $map_sth->execute;

        my $types = $map_sth->fetchrow_hashref->{Type};

        # Type will look like set('type1','type2'), so tidy it up a bit before splitting
        
        $types =~ s/set\(//;
        $types =~ s/\)//;
        $types =~ s/'//g;

        my $map;

        my $pow = 0;

        # mysql assigns the set values in consecutive powers of 2, so so shall we

        for my $type (split /,/, $types) {
            $map->{$type} = 2**$pow++;
        }

        $self->{$key} = $map;
    }
    
    return $self->{$key};
}

sub _consequences_for_set_number {
    my ($self, $set_number, $map) = @_;

    my @consequences;

    for my $term (keys %$map) {
        if ($set_number & $map->{$term}) {
            push @consequences, $OVERLAP_CONSEQUENCES{$term};
        }
    }

    return \@consequences;
}

sub _set_number_for_consequences {
    my ($self, $cons_list, $map) = @_;

    my $val = 0;

    for my $cons (@$cons_list) {
        $val |= $map->{$cons->SO_term};
    }

    return $val;
}

sub _transcript_variation_consequences_for_set_number {
    my ($self, $set_number) = @_;
    my $map = $self->_consequence_type_map('transcript_variation', 'consequence_types');
    return $self->_consequences_for_set_number($set_number, $map);
}

sub _variation_feature_consequences_for_set_number {
    my ($self, $set_number) = @_;
    my $map = $self->_consequence_type_map('variation_feature', 'consequence_type');
    return $self->_consequences_for_set_number($set_number, $map);
}

sub _transcript_variation_set_number_for_consequences {
    my ($self, $cons) = @_;
    my $map = $self->_consequence_type_map('transcript_variation', 'consequence_types');
    return $self->_set_number_for_consequences($cons, $map);
}

sub _variation_feature_set_number_for_consequences {
    my ($self, $cons) = @_;
    my $map = $self->_consequence_type_map('variation_feature', 'consequence_type');
    return $self->_set_number_for_consequences($cons, $map);
}

1;

