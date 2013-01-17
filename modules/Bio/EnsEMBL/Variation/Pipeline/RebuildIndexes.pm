=head1 LICENSE

 Copyright (c) 1999-2013 The European Bioinformatics Institute and
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
package Bio::EnsEMBL::Variation::Pipeline::RebuildIndexes;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);

sub run {

    my $self = shift;
    
    my $tables = $self->param('tables');

    my $var_dba = $self->get_species_adaptor('variation');

    my $dbc = $var_dba->dbc;

    for my $table (@$tables) {
        $dbc->do("ALTER TABLE $table ENABLE KEYS")
            or warn "Failed to enable keys on $table";
    }
}

1;

