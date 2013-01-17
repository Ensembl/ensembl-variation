=head1 LICENSE

 Copyright (c) 1999-2013 The European Bioinformatics Institute and
 Genome Research Limited.  All rights reserved.

 This software is distributed under a modified Apache license.
 For license details, please see

   http://www.ensembl.org/info/about/legal/code_licence.html

=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <dev@ensembl.org>.

 Questions may also be sent to the Ensembl help desk at
 <helpdesk@ensembl.org>.

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

