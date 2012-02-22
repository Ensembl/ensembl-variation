package Bio::EnsEMBL::Variation::Pipeline::ProteinFunction::BaseProteinFunction;

use strict;
use warnings;

use Bio::SimpleAlign;
use Bio::LocatableSeq;

use base ('Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess');

sub get_stable_id_for_md5 {
    my ($self, $md5) = @_;

    my $var_dba = $self->get_species_adaptor('variation');
    
    my $get_stable_id_sth = $var_dba->prepare(qq{
        SELECT  stable_id
        FROM    translation_mapping
        WHERE   md5 = ?
    });

    $get_stable_id_sth->execute($md5);

    my ($stable_id) = $get_stable_id_sth->fetchrow_array;

    return $stable_id;
}

sub save_predictions {
    my ($self, $pred_matrix) = @_;

    # serialize the matrix and save it to the database
    
    my $translation_md5 = $self->required_param('translation_md5');
    
    my $serialized = $pred_matrix->serialize;

    my $var_dba = $self->get_species_adaptor('variation');    

    my $dbh = $var_dba->dbc->db_handle;

    my $column = $pred_matrix->analysis;

    $column .= '_'.$pred_matrix->sub_analysis if defined $pred_matrix->sub_analysis;

    $column .= '_predictions';

    my $insert_sth = $dbh->prepare(qq{
        INSERT IGNORE INTO protein_function_predictions (translation_md5) VALUES (?)
    });

    $insert_sth->execute($translation_md5);

    my $save_sth = $dbh->prepare(qq{
        UPDATE protein_function_predictions 
        SET $column = ?
        WHERE translation_md5 = ?
    }) or die "DB error: ".$dbh->errstr;
    
    $save_sth->execute($serialized, $translation_md5);
}

1;
