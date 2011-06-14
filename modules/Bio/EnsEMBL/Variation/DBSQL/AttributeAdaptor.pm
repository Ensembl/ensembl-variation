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

1;
