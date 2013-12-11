package Bio::EnsEMBL::Variation::Pipeline::DataDumps::BaseDataDumpProcess;

use strict;

use base ('Bio::EnsEMBL::Hive::Process');


sub get_species_adaptor {
    my ($self, $species, $group) = @_;
    return $self->get_adaptor($species, $group);
}

sub get_adaptor {
    my ($self, $species, $group) = @_;
    my $dba;
    eval {
        $dba = Bio::EnsEMBL::Registry->get_DBAdaptor($species, $group);
    };
    unless (defined $dba) {
        $self->_load_registry();
        $dba = Bio::EnsEMBL::Registry->get_DBAdaptor($species, $group);
    }
    unless (defined $dba) {
        die "Failed to a get DBA for $species and group $group";
    }
    return $dba;
}

sub _load_registry {
    my ($self) = @_;
    my $reg_file = $self->param('ensembl_registry');
    Bio::EnsEMBL::Registry->load_all($reg_file, 0, 1);
    return;
}

1;

