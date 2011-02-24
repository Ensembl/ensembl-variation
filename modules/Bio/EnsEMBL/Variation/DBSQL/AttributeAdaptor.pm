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

use strict;
use warnings;

package Bio::EnsEMBL::Variation::DBSQL::AttributeAdaptor;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use Bio::EnsEMBL::Variation::OverlapConsequence;

use base qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

sub fetch_attrib_for_id {

    my ($self, $attrib_id) = @_;

    unless ($self->{attribs}) {
        
        my $attribs;

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
        }

        $self->{attribs} = $attribs;
    }

    return $self->{attribs}->{$attrib_id}->{value};
}

sub fetch_all_OverlapConsequences {
    my ($self) = @_;
    
    my @cons;
    
    for my $set (@{ $self->_fetch_sets_by_type('conseq_predicate') }) {
        push @cons, Bio::EnsEMBL::Variation::OverlapConsequence->new_fast({
            SO_accession    => $set->{SO_accession},
            SO_term         => $set->{SO_term},
            ensembl_term    => $set->{display_term},
            NCBI_term       => $set->{ncbi_term},
            feature_SO_term => $set->{feature_SO_term},
            rank            => $set->{conseq_rank},
            predicate       => $set->{conseq_predicate}
        });
    }
    
    return \@cons;
}

sub fetch_feature_type_mapping {
    my ($self) = @_;
    
    my $mapping = {};
 
    for my $set (@{ $self->_fetch_sets_by_type('ens_variant_class') }) {
        
        my $SO_term     = $set->{SO_term};
        my $class       = $set->{ens_feature_class};
        my $subtype     = $set->{ens_feature_subtype};
        my $var_type    = $set->{ens_variant_class};
        
        my $so2ens = $mapping->{$SO_term} ||= {};
            
        $so2ens->{ens_feature_class}->{$class}      = 1 if $class;
        $so2ens->{ens_feature_subtype}->{$subtype}  = 1 if $subtype;
        $so2ens->{ens_variant_class}->{$var_type}   = 1 if $var_type;
            
        my $ens2so = $mapping->{$class} ||= {};
            
        $ens2so->{SO_term}->{$SO_term} = 1 if $SO_term;
    }
    
    return $mapping;
}

sub fetch_SO_mappings {
    my ($self) = @_;
    
    my $mapping;
    
    for my $set (@{ $self->_fetch_sets_by_type('SO_term') }) {
        my $term_map = $mapping->{$set->{SO_term}} ||= {};
        $term_map->{display_term} = $set->{display_term} || $set->{SO_term};
        $term_map->{SO_accession} = $set->{SO_accession};
        
        my $acc_map = $mapping->{$set->{SO_accession}} ||= {};
        $acc_map->{SO_term} = $set->{SO_term};
        $acc_map->{display_term} = $set->{display_term} || $set->{SO_term};
    }
    
    return $mapping;
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
