=head1 LICENSE

 Copyright (c) 1999-2013 The European Bioinformatics Institute and
 Genome Research Limited.  All rights reserved.

 This software is distributed under a modified Apache license.
 For license details, please see

   http://www.ensembl.org/info/about/legal/code_licence.html

=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <dev@ensembl.org>.

 Questions may also be sent to the Ensembl help desk at
 <helpdesk@ensembl.org>.

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::AttributeAdaptor;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use Bio::EnsEMBL::Variation::OverlapConsequence;

use base qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

sub attrib_value_for_id {
    my ($self, $attrib_id) = @_;

    unless ($self->{attribs}) {
        
        my $attribs;
        my $attrib_ids;

        my $sql = qq{
            SELECT  a.attrib_id, t.code, a.value
            FROM    attrib a, attrib_type t
            WHERE   a.attrib_type_id = t.attrib_type_id
        };

        my $sth = $self->prepare($sql);

        $sth->execute;

        while (my ($attrib_id, $type, $value) = $sth->fetchrow_array) {
            $attribs->{$attrib_id}->{type}  = $type;
            $attribs->{$attrib_id}->{value} = $value;
            $attrib_ids->{$type}->{$value} = $attrib_id;
        }

        $self->{attribs}    = $attribs;
        $self->{attrib_ids} = $attrib_ids;
    }

    return defined $attrib_id ? 
        $self->{attribs}->{$attrib_id}->{value} : 
        undef;
}

sub attrib_id_for_type_value {
    my ($self, $type, $value) = @_;
    
    unless ($self->{attrib_ids}) {
        # call this method to populate the attrib hash
        $self->attrib_value_for_id;
    }
    
    return $self->{attrib_ids}->{$type}->{$value};
}

sub display_term_for_SO_term {
    my ($self, $SO_term) = @_;
    return $self->_SO_mappings->{SO_terms}->{$SO_term}->{display_term};
}

sub SO_accession_for_SO_term {
    my ($self, $SO_term) = @_;
    return $self->_SO_mappings->{SO_terms}->{$SO_term}->{SO_accession};
}

sub SO_term_for_SO_accession {
    my ($self, $SO_accession) = @_;
    return $self->_SO_mappings->{SO_accessions}->{$SO_accession}->{SO_term};
}

sub display_term_for_SO_accession {
    my ($self, $SO_accession) = @_;
    return $self->_SO_mappings->{SO_accessions}->{$SO_accession}->{display_term};
}

sub _SO_mappings {
    my ($self) = @_;
    
    unless ($self->{SO_mappings}) {
        my $mapping;
        
        for my $set (@{ $self->_fetch_sets_by_type('SO_term') }) {
            
            my $term_map = $mapping->{SO_terms}->{$set->{SO_term}} ||= {};
            $term_map->{display_term} = $set->{display_term} || $set->{SO_term};
            $term_map->{SO_accession} = $set->{SO_accession};
            
            my $acc_map = $mapping->{SO_accessions}->{$set->{SO_accession}} ||= {};
            $acc_map->{SO_term} = $set->{SO_term};
            $acc_map->{display_term} = $set->{display_term} || $set->{SO_term};
        }
        
        $self->{SO_mappings} = $mapping
    
    }
    
    return $self->{SO_mappings};
}

sub _fetch_sets_by_type {
    my ($self, $type) = @_;
    
    my $sql = qq{
        SELECT  s.attrib_set_id, t.code, a.value
        FROM    attrib a, attrib_type t, attrib_set s
        WHERE   t.attrib_type_id = a.attrib_type_id
        AND     a.attrib_id = s.attrib_id
        AND     s.attrib_set_id IN (
            SELECT s.attrib_set_id 
            FROM attrib a, attrib_type t, attrib_set s
            WHERE a.attrib_type_id = t.attrib_type_id
            AND a.attrib_id = s.attrib_id
            AND t.code = ?
        )
    };
    
    my $sth = $self->prepare($sql);
    
    $sth->execute($type);
    
    my $sets;
    
    while (my ($set_id, $type, $value) = $sth->fetchrow_array) {
        $sets->{$set_id}->{$type} = $value;
    }
    
    return [values %$sets];
}

1;
