package Bio::EnsEMBL::Variation::Pipeline::ProteinFunction::InitJobs;

use strict;

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Hive::Process;
use Digest::MD5 qw(md5_hex);

use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);

my $DEBUG = 0;

sub fetch_input {
   
    my $self = shift;
    
    my $var_dba = $self->get_species_adaptor('variation');
  
    my $fasta = $self->param('proteins_fasta');

    my $var_dbh = $var_dba->dbc->db_handle;

    # disable indexes on tables we're going to insert into as
    # this significantly speeds up inserts

    for my $table (qw(protein_info protein_position polyphen_prediction polyphen_supplementary_data sift_prediction)) {
        $var_dbh->do("ALTER TABLE $table DISABLE KEYS")
            or die "Failed to disable keys on $table: ".$var_dbh->errstr;
    }

    # create a query to add necessary new entries to the protein_position table
    
    my $add_protein_info_sth = $var_dbh->prepare_cached(qq{
        INSERT INTO protein_info (
            transcript_stable_id,
            translation_md5
        )
        VALUES (?,?)
    }) or die "DB error: ".$var_dbh->errstr;

    my $add_posn_sth = $var_dbh->prepare_cached(qq{
        INSERT INTO protein_position (
            protein_info_id, 
            position,
            amino_acid
        )
        VALUES (?,?,?)
    }) or die "DB error: ".$var_dbh->errstr;

    # loop over all toplevel slices, identifying new protein-coding transcripts

    my $tfa = $self->get_transcript_file_adaptor;

    my @transcript_stable_ids;

    if ($DEBUG) {
        @transcript_stable_ids = qw(ENST00000536460);
    }
    else {
        @transcript_stable_ids = @{ $tfa->get_all_transcript_stable_ids };
    }

    for my $transcript_stable_id (@transcript_stable_ids) {
        
        my $protein_seq = $tfa->get_protein_seq($transcript_stable_id);

        my $translation_md5 = md5_hex($protein_seq);
        
        # add an entry to the protein_info table for this transcript

        $add_protein_info_sth->execute(
            $transcript_stable_id,
            $translation_md5
        );

        my $protein_info_id = $var_dbh->last_insert_id(undef, undef, undef, undef);

        # add entries to the protein_position table for each
        # amino acid in this protein

        my @aas = split //, $protein_seq;

        my $posn = 0;

        for my $aa (@aas) {
            $posn++;
            $add_posn_sth->execute(
                $protein_info_id,
                $posn,
                $aa
            );
        }
    }

    # we're done inserting into protein_info and protein_position now, so enable the indexes again

    for my $table (qw(protein_info protein_position)) {
        $var_dbh->do("ALTER TABLE $table ENABLE KEYS")
            or die "Failed to enable keys on $table: ".$var_dbh->errstr;
    }

    my @output_ids = map { {transcript_stable_id => $_} } @transcript_stable_ids;

    $self->param('polyphen_output_ids', \@output_ids);
    $self->param('sift_output_ids', \@output_ids);
}

sub write_output {
    my $self = shift;
    
    $self->dataflow_output_id({}, 1); # rebuild polyphen indexes
    $self->dataflow_output_id($self->param('polyphen_output_ids'), 2);
    $self->dataflow_output_id({}, 3); # rebuild sift indexes
    $self->dataflow_output_id($self->param('sift_output_ids'), 4);
}

1;
