package Bio::EnsEMBL::Variation::Pipeline::ProteinFunction::InitJobs;

use strict;

use Bio::EnsEMBL::Registry;

use Bio::EnsEMBL::Variation::Pipeline::ProteinFunction::Constants qw(FULL UPDATE NONE);

use Digest::MD5 qw(md5_hex);

use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);

my $DEBUG = 0;

sub fetch_input {
   
    my $self = shift;
    
    my $var_dba  = $self->get_species_adaptor('variation');
    my $core_dba = $self->get_species_adaptor('core');
  
    my $fasta = $self->required_param('proteins_fasta');

    my $sift_run_type       = $self->required_param('sift_run_type');
    my $polyphen_run_type   = $self->required_param('polyphen_run_type');
    my $include_lrg         = $self->param('include_lrg');
    my $use_existing_table  = $self->param('use_existing_table');

    my $var_dbh = $var_dba->dbc->db_handle;
   
    my $get_existing_sth = $var_dbh->prepare(qq{
        SELECT  translation_md5, transcript_stable_id, polyphen_predictions, sift_predictions
        FROM    protein_function_predictions
        WHERE   translation_stable_id = ?
    });

    my $delete_preds_sth = $var_dbh->prepare(qq{
        DELETE FROM protein_function_predictions
        WHERE   translation_stable_id = ?
    });

    my $add_preds_sth = $var_dbh->prepare_cached(qq{
        INSERT INTO protein_function_predictions (
            translation_stable_id,
            transcript_stable_id,
            translation_md5
        )
        VALUES (?,?,?)
    }) or die "DB error: ".$var_dbh->errstr;

    my $fix_transcript_stable_id = $var_dbh->prepare(qq{
        UPDATE  protein_function_predictions
        SET     transcript_stable_id = ?
        WHERE   translation_stable_id = ?
    });

    my $ta = $core_dba->get_TranscriptAdaptor or die "Failed to get transcript adaptor";
    my $ga = $core_dba->get_GeneAdaptor or die "Failed to get gene adaptor";
    my $sa = $core_dba->get_SliceAdaptor or die "Failed to get slice adaptor";

    # fetch all the regular genes

    my $transcripts;

    if ($DEBUG) {
        $transcripts = $ga->fetch_all_by_external_name('BRCA2')->[0]->get_all_Transcripts;
    }
    else {
        for my $slice (@{ $sa->fetch_all('toplevel', undef, 1, undef, ($include_lrg ? 1 : undef)) }) {
            for my $gene (@{ $slice->get_all_Genes(undef, undef, 1) }) {
                for my $transcript (@{ $gene->get_all_Transcripts }) {
                    if (my $translation = $transcript->translation) {
                        push @$transcripts, $transcript;
                    }
                }
            }
        }
    }

    # get a transcript file adaptor over these transcripts

    my $tfa;

    if ($use_existing_table) {
        $tfa = $self->get_transcript_file_adaptor;
    }
    else {
        # create a new fasta dump of all the transcripts
        $tfa = $self->get_transcript_file_adaptor($transcripts);
    }

    my @sift_stable_ids;
    my @polyphen_stable_ids;

    for my $transcript (@$transcripts) {
        
        my $translation_stable_id = $transcript->translation->stable_id;

        my $translation_md5;
        
        unless ($use_existing_table) {
            
            # check the table to see if we need to run anything

            $get_existing_sth->execute($translation_stable_id);

            my (
                $existing_md5, 
                $existing_transcript_stable_id, 
                $pph_preds, 
                $sift_preds
            ) = $get_existing_sth->fetchrow_array;

            #$translation_md5 = $tfa->get_translation_md5($translation_stable_id);
            $translation_md5 = md5_hex($transcript->translation->seq);

            if ($existing_md5) {

                if ($existing_md5 eq $translation_md5) {

                    # the protein hasn't changed, so we'll only run the analyses for
                    # the protein if we're doing a full new run, or if we don't have any
                    # existing predictions

                    if ($sift_run_type == FULL || length($sift_preds) < 1) {
                        push @sift_stable_ids, $translation_stable_id;
                    }

                    if ($polyphen_run_type == FULL || length($pph_preds) < 1) {
                        push @polyphen_stable_ids, $translation_stable_id;      
                    }
                    
                    # check that the transcript_stable_id for this translation is correct

                    unless ($existing_transcript_stable_id eq $transcript->stable_id) {
                        $fix_transcript_stable_id->execute($transcript->stable_id, $translation_stable_id);
                    }

                    # next out of the loop here as we'll keep the same entry in
                    # the protein_function_predictions table

                    next;
                }
                else {
                    # the translation has changed, so delete the existing entry
                    
                    $delete_preds_sth->execute($translation_stable_id);
                }
            }
        }

        # this is a new or altered protein so we need to run both tools, unless their run type is NONE
        
        push @sift_stable_ids, $translation_stable_id unless $sift_run_type == NONE;
        push @polyphen_stable_ids, $translation_stable_id unless $polyphen_run_type == NONE;      

        # and we need to add this protein to the protein_function_predictions table

        unless ($use_existing_table) {
            $add_preds_sth->execute(
                $translation_stable_id,
                $transcript->stable_id,
                $translation_md5
            );
        }
    }
    
    # and set up the necessary sift and polyphen jobs

    my @polyphen_output_ids = map { {translation_stable_id => $_} } @polyphen_stable_ids;
    my @sift_output_ids = map { {translation_stable_id => $_} } @sift_stable_ids;

    $self->param('polyphen_output_ids', \@polyphen_output_ids);
    $self->param('sift_output_ids', \@sift_output_ids);
}

sub write_output {
    my $self = shift;
    
    unless ($self->param('polyphen_run_type') == NONE) {
        $self->dataflow_output_id($self->param('polyphen_output_ids'), 2);
    }

    unless ($self->param('sift_run_type') == NONE) {
        $self->dataflow_output_id($self->param('sift_output_ids'), 3);
    }
}

1;
