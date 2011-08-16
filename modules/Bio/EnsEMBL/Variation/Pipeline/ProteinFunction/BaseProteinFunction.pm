package Bio::EnsEMBL::Variation::Pipeline::ProteinFunction::BaseProteinFunction;

use strict;
use warnings;

use Bio::SimpleAlign;
use Bio::LocatableSeq;

use Bio::EnsEMBL::Variation::Utils::ProteinFunctionUtils qw(
    @ALL_AAS
    $NO_PREDICTION
    prediction_to_short 
    compress_prediction_string
);

use base ('Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess');

sub save_predictions {
    my ($self, $tool, $preds) = @_;

    my $translation_stable_id = $self->required_param('translation_stable_id');

    my $peptide = $self->get_transcript_file_adaptor->get_translation_seq($translation_stable_id);

    # build up our prediction string, making sure we 
    # have an entry for every possible position and amino
    # acid, even if it's just the dummy $NO_PREDICTION value

    my $pred_str;

    for my $pos (1 .. length($peptide)) {
        
        for my $aa (@ALL_AAS) {
            
            my $pred = undef;

            my $prediction  = $preds->{$pos}->{$aa}->{prediction};
            my $score       = $preds->{$pos}->{$aa}->{score};
            
            if (defined $prediction && defined $score) {
                $pred = prediction_to_short($tool, $prediction, $score);
            }

            $pred_str .= defined $pred ? $pred : $NO_PREDICTION;
        }
    }

    # compress the string

    my $gzipped = compress_prediction_string($pred_str);

    # and save it to the database

    my $var_dba = $self->get_species_adaptor('variation');    

    my $dbh = $var_dba->dbc->db_handle;
   
    my $col = $tool.'_predictions';

    my $save_sth = $dbh->prepare(qq{
        UPDATE protein_function_predictions SET $col = ?
        WHERE translation_stable_id = ?
    }) or die "DB error: ".$dbh->errstr;
    
    $save_sth->execute($gzipped, $translation_stable_id);
}

1;
