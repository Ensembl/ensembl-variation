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
package Bio::EnsEMBL::Variation::Pipeline::FinishVariationClass;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);

my $DEBUG = 0;

sub run {
    my $self = shift;
    
    my $temp_var_table = $self->required_param('temp_var_table');
    
    my $temp_var_feat_table = $self->required_param('temp_var_feat_table');
  
    my $var_dba = $self->get_species_adaptor('variation');
    
    my $dbc = $var_dba->dbc;
    
    # copy the attribs across to the real tables

    $dbc->do(qq{
        UPDATE  variation_feature vf, $temp_var_feat_table tvf
        SET     vf.class_attrib_id = tvf.class_attrib_id
        WHERE   vf.variation_feature_id = tvf.variation_feature_id
    });

    $dbc->do(qq{
        UPDATE  variation v, $temp_var_table tv
        SET     v.class_attrib_id = tv.class_attrib_id
        WHERE   v.variation_id = tv.variation_id
    });
    
    # and get rid of the temp tables
    $dbc->do(qq{DROP TABLE $temp_var_table});
    $dbc->do(qq{DROP TABLE $temp_var_feat_table});
    
    return;
}

1;

