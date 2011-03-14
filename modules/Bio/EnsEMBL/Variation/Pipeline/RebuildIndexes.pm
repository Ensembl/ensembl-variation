package Bio::EnsEMBL::Variation::Pipeline::RebuildIndexes;

use strict;

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Hive::Process;

use base ('Bio::EnsEMBL::Hive::Process');

sub run {

    my $self = shift;

    my $reg_file = $self->param('ensembl_registry')
        or die "ensembl_registry is a required parameter";
    
    my $species = $self->param('species')
        or die "species is a required parameter";
    
    my $tables = $self->param('tables')
        or die "tables is a required parameter";

    my $reg = 'Bio::EnsEMBL::Registry';
    
    $reg->load_all($reg_file, 0, 1);
  
    my $var_dba = $reg->get_DBAdaptor($species, 'variation')
        or die "failed to get variation DBA for $species";

    my $dbh = $var_dba->dbc->db_handle;

    for my $table (@$tables) {
        $dbh->do("ALTER TABLE $table ENABLE KEYS")
            or warn "Failed to enable keys on $table: ".$dbh->errstr;
    }
}

1;

