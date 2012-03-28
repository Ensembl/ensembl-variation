package Bio::EnsEMBL::Variation::Pipeline::ProteinFunction::InitJobs;

use strict;

use Bio::EnsEMBL::Registry;

use Bio::EnsEMBL::Variation::Pipeline::ProteinFunction::Constants qw(FULL UPDATE NONE);

use Digest::MD5 qw(md5_hex);

use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);

sub fetch_input {
   
    my $self = shift;
    
    my $sift_run_type   = $self->required_param('sift_run_type');
    my $pph_run_type    = $self->required_param('pph_run_type');
    my $include_lrg     = $self->param('include_lrg');
    
    my $core_dba = $self->get_species_adaptor('core');
    my $var_dba  = $self->get_species_adaptor('variation');

    # fetch all the transcripts from the core DB

    my @transcripts;

    if ($self->param('debug_mode')) {
        my $ga = $core_dba->get_GeneAdaptor or die "Failed to get gene adaptor";
        
        @transcripts = grep { $_->translation } @{ $ga->fetch_all_by_external_name('BRCA1')->[0]->get_all_Transcripts };
    }
    else {
        my $sa = $core_dba->get_SliceAdaptor or die "Failed to get slice adaptor";
        
        for my $slice (@{ $sa->fetch_all('toplevel', undef, 1, undef, ($include_lrg ? 1 : undef)) }) {
            for my $gene (@{ $slice->get_all_Genes(undef, undef, 1) }) {
                for my $transcript (@{ $gene->get_all_Transcripts }) {
                    if (my $translation = $transcript->translation) {
                        push @transcripts, $transcript;
                    }
                }
            }
        }
    }

    # store a table mapping each translation stable ID to its corresponding MD5

    $var_dba->dbc->do(qq{DROP TABLE IF EXISTS translation_mapping});

    $var_dba->dbc->do(qq{
        CREATE TABLE translation_mapping (
            stable_id   VARCHAR(255),
            md5         CHAR(32),
            PRIMARY KEY (stable_id),
            KEY md5_idx (md5)
        )
    });

    my $add_mapping_sth = $var_dba->prepare(qq{
        INSERT IGNORE INTO translation_mapping (stable_id, md5) VALUES (?,?)
    });

    my $add_md5_sth = $var_dba->prepare(qq{
        INSERT IGNORE INTO translation_md5 (translation_md5) VALUES (?)
    });

    # build a hash mapping MD5s for our set of translations to their peptide sequences
    # and also write each stable ID - MD5 mapping to the database

    my %unique_translations;

    for my $tran (@transcripts) {
        
        my $tl = $tran->translation;

        my $seq = $tl->seq;
        
        my $md5 = md5_hex($seq);
        
        $unique_translations{$md5} = $seq;

        $add_mapping_sth->execute($tl->stable_id, $md5);
    }

    # work out which translations we need to run for which analysis
    
    my @translation_md5s = keys %unique_translations;
    
    my @sift_md5s;
    my @pph_md5s;

    # if we're doing full runs, then we just use all the translations

    @sift_md5s = @translation_md5s if $sift_run_type == FULL;
    @pph_md5s  = @translation_md5s if $pph_run_type == FULL;

    # if we're updating we need to check which translations already have predictions
        
    if ($sift_run_type == UPDATE || $pph_run_type == UPDATE) {
        
        my $var_dbh = $var_dba->dbc->db_handle;
    
        my $get_existing_sth = $var_dbh->prepare(qq{
            SELECT  t.translation_md5, a.value, p.prediction_matrix
            FROM    translation_md5 t, attrib a, protein_function_predictions p
            WHERE   t.translation_md5_id = p.translation_md5_id 
            AND     a.attrib_id = p.analysis_attrib_id
        });

        # store the set of existing MD5s in a hash

        $get_existing_sth->execute;

        my $existing_md5s;

        while ( my ($md5, $analysis, $preds) = $get_existing_sth->fetchrow_array ) {
            # there are 2 polyphen analyses, but we only want to track one
            $analysis = 'pph' if $analysis =~ /polyphen/;
            # just record true if we already have predictions for each translation
            $existing_md5s->{$md5}->{$analysis} = length($preds) > 0;   
        }

        # exclude any translations MD5s we find in the database

        @sift_md5s = grep { ! $existing_md5s->{$_}->{sift} } @translation_md5s if $sift_run_type == UPDATE;
        @pph_md5s  = grep { ! $existing_md5s->{$_}->{pph} } @translation_md5s if $pph_run_type == UPDATE;
    }
    
    # create a FASTA dump of all the necessary translation sequences

    # work out the set of all unique MD5s for sift and polyphen

    my %required_md5s = map { $_ => 1 } (@sift_md5s, @pph_md5s);

    my $fasta = $self->required_param('fasta_file');

    open my $FASTA, ">$fasta" or die "Failed to open $fasta for writing";

    # get rid of any existing index file

    if (-e "$fasta.fai") {
        unlink "$fasta.fai" or die "Failed to delete fasta index file";
    }
    
    # dump out each peptide sequence

    for my $md5 (keys %required_md5s) {
        my $seq = $unique_translations{$md5};
        $seq =~ s/(.{80})/$1\n/g;
        # get rid of any trailing newline
        chomp $seq;
        print $FASTA ">$md5\n$seq\n";
    }

    close $FASTA;

    # set up our list of output ids

    $self->param('pph_output_ids',  [ map { {translation_md5 => $_} } @pph_md5s ]);
    $self->param('sift_output_ids', [ map { {translation_md5 => $_} } @sift_md5s ]);
}

sub write_output {
    my $self = shift;
    
    unless ($self->param('pph_run_type') == NONE) {
        $self->dataflow_output_id($self->param('pph_output_ids'), 2);
    }

    unless ($self->param('sift_run_type') == NONE) {
        $self->dataflow_output_id($self->param('sift_output_ids'), 3);
    }
}

1;
