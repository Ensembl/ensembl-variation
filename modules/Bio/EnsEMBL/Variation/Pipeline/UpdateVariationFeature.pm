package Bio::EnsEMBL::Variation::Pipeline::UpdateVariationFeature;

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
    
    my $reg = 'Bio::EnsEMBL::Registry';
    
    $reg->load_all($reg_file, 0, 1);
  
    my $var_dba = $reg->get_DBAdaptor($species, 'variation')
        or die "failed to get variation DBA for $species";

    my $dbh = $var_dba->dbc->db_handle;
    
    # first set the default consequence type

    $dbh->do(qq{
        UPDATE  variation_feature
        SET     consequence_type = 'intergenic_variant'
    }) or die "Failed to reset consequence_type on variation_feature";

    # create a temp table (dropping it if it exists)

    $dbh->do(qq{DROP TABLE IF EXISTS variation_feature_with_tv})
        or die "Failed to drop pre-existing temp table";
    
    $dbh->do(qq{CREATE TABLE variation_feature_with_tv LIKE variation_feature})
        or die "Failed to create temp table";

    # concatenate the consequence types from transcript_variation 

    $dbh->do(qq{
        INSERT INTO variation_feature_with_tv (variation_feature_id, consequence_type)
        SELECT  variation_feature_id, GROUP_CONCAT(consequence_types) 
        FROM    transcript_variation 
        GROUP BY variation_feature_id
    }) or die "Populating temp table failed";

    # update variation_feature
    
    $dbh->do(qq{
        UPDATE  variation_feature vf, variation_feature_with_tv vftv
        SET     vf.consequence_type = vftv.consequence_type
        WHERE   vf.variation_feature_id = vftv.variation_feature_id
    }) or die "Failed to update vf table";

    # and get rid of our temp table

    $dbh->do(qq{DROP TABLE variation_feature_with_tv})
        or die "Failed to drop temp table";

}


1;

