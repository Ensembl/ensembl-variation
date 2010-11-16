package Bio::EnsEMBL::Variation::Pipeline::RebuildConsequenceIndexes;

use strict;

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Hive::Process;

use base ('Bio::EnsEMBL::Hive::Process');

sub run {
    
    my $reg_file = $self->param('ensembl_registry')
        or die "ensembl_registry is a required parameter";
    
    my $species = $self->param('species')
        or die "species is a required parameter";
 
    my $reg = 'Bio::EnsEMBL::Registry';
    
    $reg->load_all($reg_file, 0, 1);
  
    my $var_dba = $reg->get_DBAdaptor($species, 'variation')
        or die "failed to get variation DBA for $species";

    my $dbh = $var_dba->dbc->db_handle;

    $dbh->do(q{
        CREATE INDEX feature_type_idx       ON vf_overlap(feature_type_id);
        CREATE INDEX variation_feature_idx  ON vf_overlap(variation_feature_id);
        CREATE INDEX feature_stable_idx     ON vf_overlap(feature_stable_id);
        CREATE INDEX vf_overlap_idx         ON vf_overlap_allele(vf_overlap_id);
        CREATE INDEX type_val_idx           ON vf_overlap_attrib(attrib_type_id, value(40));
        CREATE INDEX val_only_idx           ON vf_overlap_attrib(value(40));
        CREATE INDEX overlap_idx            ON vf_overlap_attrib(vf_overlap_id);
        CREATE INDEX type_val_idx           ON vf_overlap_allele_attrib(attrib_type, value(40));
        CREATE INDEX val_only_idx           ON vf_overlap_allele_attrib(value(40));
        CREATE INDEX allele_idx             ON vf_overlap_allele_attrib(vf_overlap_allele_id);
    }) or die "Failed to rebuild indexes: ".$dbh->errstr;

}

1;
