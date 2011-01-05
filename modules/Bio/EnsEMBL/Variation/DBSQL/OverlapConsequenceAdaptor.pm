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

package Bio::EnsEMBL::Variation::DBSQL::OverlapConsequenceAdaptor;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use Bio::EnsEMBL::Variation::OverlapConsequence;

use base qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

sub _tables {
    return (
        ['overlap_consequence'],
    );
}

sub _columns {
    return qw (
        overlap_consequence_id
        SO_term
        SO_id
        ensembl_term
        NCBI_term
        feature_SO_term
        rank
        predicate
    );
}

sub _objs_from_sth {
    my ($self, $sth) = @_;
    
    my @objs;
    
    while (my $row = $sth->fetchrow_hashref) {
        
        # rename overlap_consequence_id to dbID
        $row->{dbID} = $row->{overlap_consequence_id};
        delete $row->{overlap_consequence_id};
        
        push @objs, Bio::EnsEMBL::Variation::OverlapConsequence->new_fast($row);
    }
    
    return \@objs;
}

1;
