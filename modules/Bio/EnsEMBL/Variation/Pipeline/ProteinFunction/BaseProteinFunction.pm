package Bio::EnsEMBL::Variation::Pipeline::ProteinFunction::BaseProteinFunction;

use strict;
use warnings;

use Bio::SimpleAlign;
use Bio::LocatableSeq;

use base ('Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess');

sub save_predictions {
    my ($self, $pred_matrix) = @_;

    # serialize the matrix and save it to the database
    
    my $translation_stable_id = $self->required_param('translation_stable_id');
    
    my $binary_string = $pred_matrix->serialize;

    my $var_dba = $self->get_species_adaptor('variation');    

    my $dbh = $var_dba->dbc->db_handle;
   
    my $col = $pred_matrix->analysis.'_predictions';

    my $save_sth = $dbh->prepare(qq{
        UPDATE protein_function_predictions SET $col = ?
        WHERE translation_stable_id = ?
    }) or die "DB error: ".$dbh->errstr;
    
    $save_sth->execute($binary_string, $translation_stable_id);
}

1;
