package Bio::EnsEMBL::Variation::Pipeline::ProteinFunction::RunSift;

use strict;

use Digest::MD5 qw(md5_hex);
use File::Path qw(make_path remove_tree);
use Data::Dumper;

use Bio::EnsEMBL::Registry;

use base ('Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess');

my $DEBUG   = 0;

my $MEDIAN_CUTOFF = 2.75; # as per README

sub run {
    my $self = shift;

    my $transcript_stable_id    = $self->param('transcript_stable_id');
    my $fasta                   = $self->param('proteins_fasta');
    my $sift_dir                = $self->param('sift_dir');
    my $ncbi_dir                = $self->param('ncbi_dir');
    my $blastdb                 = $self->param('blastdb');

    my $md5 = md5_hex($transcript_stable_id);
    my $dir = substr($md5, 0, 2);
    my $output_dir = "$sift_dir/working/$dir/$transcript_stable_id";
    my $tarball = 'scratch.tgz';

    unless (-d $output_dir) {
        my $err;
        make_path($output_dir, {error => \$err});
        die "make_path failed: ".Dumper($err) if $err && @$err;
    }

    chdir $output_dir or die "Failed to chdir to $output_dir";

    #my $root_file  = "$output_dir/$transcript_stable_id";
    my $root_file  = "$transcript_stable_id";
    my $fasta_file = "$root_file.fa";
    my $aln_file   = "$root_file.alignedfasta";
    my $res_file   = "$root_file.SIFTprediction";
    my $subs_file  = "$root_file.subst";

    if (-e "$output_dir/$tarball") {
        system("tar zxvf $tarball > /dev/null") == 0
            or die "Failed to untar $output_dir/$tarball: $!";
    }
    
    # set necessary environment variables for sift

    $ENV{NCBI}       = $ncbi_dir;
    $ENV{BLIMPS_DIR} = $sift_dir.'/blimps';
    $ENV{SIFT_DIR}   = $sift_dir;   
    $ENV{tmpdir}     = $output_dir;   
    
    # fetch our proteins 

    my $peptide = $self->get_transcript_file_adaptor->get_protein_seq($transcript_stable_id);

    my $alignment_ok = 1;

    $self->dbc->disconnect_when_inactive(1);

    unless (-e $aln_file) {
        
        # we need to do a multiple alignment

        # first create a fasta file for the protein of this transcript (there should always 
        # be one because we only create a job for transcripts with a translation)
    
        open (FASTA_FILE, ">$fasta_file");

        my $pep_str = $peptide;

        $pep_str =~ s/(.{80})/$1\n/g;

        print FASTA_FILE '>'.$transcript_stable_id,"\n$pep_str\n";

        close FASTA_FILE;

        # and run the alignment program

        my $cmd = "$sift_dir/bin/ensembl_seqs_chosen_via_median_info.csh $fasta_file $blastdb $MEDIAN_CUTOFF";

        my $exit_code = system($cmd);
        
        if ($exit_code == 0) {
            $alignment_ok = 1;
        }
        else {
            # the alignment failed for some reason, what to do?
            warn "Alignment for $transcript_stable_id failed - cmd: $cmd";
            $alignment_ok = 0;
        }
    }

    if ($alignment_ok) {

        # work out the sift score for each possible amino acid substitution

        unless (-e $subs_file) {

            # create our substitution file

            my @aas = split //, $peptide;

            my @all_aas = qw(A C D E F G H I K L M N P Q R S T V W Y);

            my $pos = 0;

            open SUBS, ">$subs_file" or die "Failed to open $subs_file: $!";

            for my $ref (@aas) {
                $pos++;
                for my $alt (@all_aas) {
                    unless ($ref eq $alt) {
                        print SUBS $ref.$pos.$alt."\n";
                    }
                }
            }

            close SUBS;
        }

        # and run sift on it

        my $cmd = "$sift_dir/bin/info_on_seqs $aln_file $subs_file $res_file";

        system($cmd) == 0 or die "Failed to run $cmd: $?";

        $self->dbc->disconnect_when_inactive(0);

        # parse and store the results 
        
        my $var_dba = $self->get_species_adaptor('variation');
  
        my $dbh = $var_dba->dbc->db_handle;

        # define the necessary SQL statements

        my $get_pos_sth = $dbh->prepare_cached(qq{
            SELECT  pp.protein_position_id, pp.position
            FROM    protein_position pp, protein_info pi
            WHERE   pp.protein_info_id = pi.protein_info_id
            AND     pi.transcript_stable_id = ?
        }) or die "DB error: ".$dbh->errstr;

        my $add_sift_pos_info_sth = $dbh->prepare_cached(qq{
            UPDATE  protein_position
            SET     sift_median_conservation = ?, sift_num_sequences_represented = ?
            WHERE   protein_position_id = ?
        }) or die "DB error: ".$dbh->errstr;

        my $save_sth = $dbh->prepare_cached(qq{
            INSERT INTO sift_prediction (
                protein_position_id,
                amino_acid,
                score
            ) 
            VALUES (?,?,?)
        }) or die "DB error: ".$dbh->errstr;

        open (RESULTS, "<$res_file") or die "Failed to open $res_file: $!";

        # fetch the protein_position_ids for this protein

        $get_pos_sth->execute($transcript_stable_id);

        my %pos_ids;

        while (my ($pos_id, $pos) = $get_pos_sth->fetchrow_array) {
            $pos_ids{$pos} = $pos_id;
        }

        # parse the results file

        my %pos_has_sift_info;

        while (<RESULTS>) {

            chomp;

            next if /WARNING/;
            next if /NOT SCORED/;

            my ($subst, $prediction, $score, $median_cons, $num_seqs, $blocks) = split;

            my ($ref_aa, $pos, $alt_aa) = $subst =~ /([A-Z])(\d+)([A-Z])/;

            next unless $ref_aa && $alt_aa && defined $pos;

            my $pos_id = $pos_ids{$pos} || die "No protein_position entry for $transcript_stable_id pos $pos";

            unless ($pos_has_sift_info{$pos_id}) {
                # add the sift position specific information

                $add_sift_pos_info_sth->execute(
                    $median_cons,
                    $num_seqs,
                    $pos_id
                );

                $pos_has_sift_info{$pos_id}++;
            }

            # save the prediction for this specific possible substitution

            $save_sth->execute(
                $pos_id,
                $alt_aa,
                $score            
            );
        }
    }

    $self->dbc->disconnect_when_inactive(0);

    # tar up the files

    system("tar --remove-files --exclude *.tgz -czvf $tarball * > /dev/null") == 0
        or die "Failed to create $output_dir/$tarball: $!";
}
