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

