=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut


=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

 Questions may also be sent to the Ensembl help desk at
 <http://www.ensembl.org/Help/Contact>.

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

sub attrib_id_for_type_code {
    my ($self, $type) = @_;

    unless ($self->{attrib_types}) {
        
        my $attrib_types;

        my $sql = qq{
            SELECT  t.attrib_type_id, t.code, t.name, t.description
            FROM    attrib_type t
        };

        my $sth = $self->prepare($sql);

        $sth->execute;

        while (my ($attrib_type_id, $code, $name, $description ) = $sth->fetchrow_array) {
            $attrib_types->{$code}->{attrib_type_id} = $attrib_type_id;
            $attrib_types->{$code}->{name}           = ($name eq '') ? $code : $name;
            $attrib_types->{$code}->{description}    = $description;
        }

        $self->{attrib_types}  = $attrib_types;
    }

    return defined $type ? 
        $self->{attrib_types}->{$type}->{attrib_type_id} :
        undef;
}

sub attrib_type_name_for_attrib_type_code {
    my ($self, $type) = @_;

    unless ($self->{attrib_types}) {
        # call this method to populate the attrib_types hash
        $self->attrib_id_for_type_code;
    }

    return defined $type ?
        $self->{attrib_types}->{$type}->{name} :
        undef;
}

sub attrib_type_code_for_attrib_type_id {
    my ($self, $a_type_id) = @_;

    unless ($self->{attrib_types_by_ids}) {

        my $attrib_types;

        my $sql = qq{
            SELECT  t.attrib_type_id, t.code
            FROM    attrib_type t
        };

        my $sth = $self->prepare($sql);

        $sth->execute;

        while (my ($attrib_type_id, $code ) = $sth->fetchrow_array) {
            $attrib_types->{$attrib_type_id}  = $code ;
        }

        $self->{attrib_types_by_ids}  = $attrib_types;
    }

    return defined $a_type_id ?
        $self->{attrib_types_by_ids}->{$a_type_id} :
        undef;
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
