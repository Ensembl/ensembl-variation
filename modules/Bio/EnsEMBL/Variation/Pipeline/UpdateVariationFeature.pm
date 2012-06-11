package Bio::EnsEMBL::Variation::Pipeline::UpdateVariationFeature;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);

sub run {

    my $self = shift;

    my $var_dba = $self->get_species_adaptor('variation');
    
    my $dbc = $var_dba->dbc;
    
    # first set the default consequence type

    $dbc->do(qq{
        UPDATE  variation_feature
        SET     consequence_types = 'intergenic_variant'
    }) or die "Failed to reset consequence_types on variation_feature";

    # create a temp table (dropping it if it exists)

    my $temp_table = 'variation_feature_with_tv';

    $dbc->do(qq{DROP TABLE IF EXISTS $temp_table})
        or die "Failed to drop pre-existing temp table";
    
    $dbc->do(qq{CREATE TABLE $temp_table LIKE variation_feature})
        or die "Failed to create temp table";

    # concatenate the consequence types from transcript_variation 

    $dbc->do(qq{
        INSERT INTO $temp_table (variation_feature_id, consequence_types)
        SELECT  variation_feature_id, GROUP_CONCAT(DISTINCT(consequence_types)) 
        FROM    transcript_variation 
        GROUP BY variation_feature_id
    }) or die "Populating temp table failed";

    # update variation_feature
    
    $dbc->do(qq{
        UPDATE  variation_feature vf, $temp_table tvf
        SET     vf.consequence_types = tvf.consequence_types
        WHERE   vf.variation_feature_id = tvf.variation_feature_id
    }) or die "Failed to update vf table";

    # and get rid of our temp table

    $dbc->do(qq{DROP TABLE $temp_table})
        or die "Failed to drop temp table";

}

1;

