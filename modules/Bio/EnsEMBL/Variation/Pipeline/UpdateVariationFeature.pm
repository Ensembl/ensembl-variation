=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut


=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

 Questions may also be sent to the Ensembl help desk at
 <http://www.ensembl.org/Help/Contact>.

=cut
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

    ## remove unneccessary non-null columns (for EGenomes)

    $dbc->do(qq{ALTER TABLE $temp_table  drop seq_region_id,       drop variation_id ,
                                         drop seq_region_start,    drop seq_region_end,
                                         drop seq_region_strand,   drop source_id,
                                         drop map_weight })
        or die "Failed to alter temp table";

   
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

