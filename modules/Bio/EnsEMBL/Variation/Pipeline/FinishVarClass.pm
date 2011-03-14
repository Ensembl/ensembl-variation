package Bio::EnsEMBL::Variation::Pipeline::FinishVarClass;

use strict;

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Hive::Process;

use base ('Bio::EnsEMBL::Hive::Process');

my $DEBUG = 0;

sub run {
    
    my $self = shift;

    my $reg_file = $self->param('ensembl_registry')
        or die "ensembl_registry is a required parameter";

    my $species = $self->param('species')
        or die "species is a required parameter";
    
    my $temp_var_table = $self->param('temp_var_table')
        or die "temp_var_table is required";
    
    my $temp_var_feat_table = $self->param('temp_var_feat_table')
        or die "temp_var_feat_table is required";

    my $reg = 'Bio::EnsEMBL::Registry';
    
    $reg->load_all($reg_file, 0, 1);
  
    my $var_dba = $reg->get_DBAdaptor($species, 'variation')
        or die "failed to get variation DBA for $species";

    my $dbh = $var_dba->dbc->db_handle;

    # copy the attribs across to the real tables

    $dbh->do(qq{
        UPDATE  variation_feature vf, $temp_var_feat_table tvf
        SET     vf.class_attrib_id = tvf.class_attrib_id
        WHERE   vf.variation_feature_id = tvf.variation_feature_id
    });

    $dbh->do(qq{
        UPDATE  variation v, $temp_var_table tv
        SET     v.class_attrib_id = tv.class_attrib_id
        WHERE   v.variation_id = tv.variation_id
    });
    
    # and get rid of the temp tables

    $dbh->do(qq{DROP TABLE $temp_var_table});
    $dbh->do(qq{DROP TABLE $temp_var_feat_table});
}

1;

