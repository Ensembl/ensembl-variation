package Bio::EnsEMBL::Variation::Pipeline::ProteinFunction::RunSift;

use strict;

use Digest::MD5 qw(md5_hex);
use File::Path qw(make_path remove_tree);
use Data::Dumper;

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix qw(@ALL_AAS);
use Bio::EnsEMBL::Variation::Utils::ComparaUtils qw(dump_alignment_for_sift);

use base ('Bio::EnsEMBL::Variation::Pipeline::ProteinFunction::BaseProteinFunction');

my $DEBUG   = 0;

my $MEDIAN_CUTOFF = 2.75; # as per README

sub run {
    my $self = shift;

    my $translation_stable_id   = $self->required_param('translation_stable_id');
    my $sift_dir                = $self->required_param('sift_dir');
    my $ncbi_dir                = $self->required_param('ncbi_dir');
    my $blastdb                 = $self->required_param('blastdb');

    my $md5 = md5_hex($translation_stable_id);
    my $dir = substr($md5, 0, 2);
    my $output_dir = "$sift_dir/working/$dir/$translation_stable_id";
    my $tarball = 'scratch.tgz';

    unless (-d $output_dir) {
        my $err;
        make_path($output_dir, {error => \$err});
        die "make_path failed: ".Dumper($err) if $err && @$err;
    }

    chdir $output_dir or die "Failed to chdir to $output_dir";

    my $root_file  = "$translation_stable_id";
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
    
    # fetch our protein 

    my $tfa = $self->get_transcript_file_adaptor;

    my $peptide = $tfa->get_translation_seq($translation_stable_id);

    die "No protein sequence for $translation_stable_id" unless $peptide && (length($peptide) > 0);

    my $alignment_ok = 1;

    unless (-e $aln_file) {
        
        # we need to get the multiple alignment
        
        if ($self->param('use_compara')) {
            
            eval {
                dump_alignment_for_sift($translation_stable_id, $aln_file);
            };

            if ($@) {
                warn "Failed to get a compara alignment for $translation_stable_id: $@";
                $alignment_ok = 0;
            }
        }
        else {

            # do the alignment ourselves
            
            # first create a fasta file for the protein sequence

            open (FASTA_FILE, ">$fasta_file");

            print FASTA_FILE $tfa->get_translation_fasta($translation_stable_id);

            close FASTA_FILE;

            # and run the alignment program

            $self->dbc->disconnect_when_inactive(1);

            my $cmd = "$sift_dir/bin/ensembl_seqs_chosen_via_median_info.csh $fasta_file $blastdb $MEDIAN_CUTOFF";

            $self->dbc->disconnect_when_inactive(0);

            my $exit_code = system($cmd);
            
            if ($exit_code == 0) {
                $alignment_ok = 1;
            }
            else {
                # the alignment failed for some reason, what to do?
                warn "Alignment for $translation_stable_id failed - cmd: $cmd";
                $alignment_ok = 0;
            }
        }
    }

    if ($alignment_ok) {

        # work out the sift score for each possible amino acid substitution

        unless (-e $subs_file) {

            # create our substitution file

            my $pos = 0;

            open SUBS, ">$subs_file" or die "Failed to open $subs_file: $!";
            
            my @aas = split //, $peptide;

            for my $ref (@aas) {
                $pos++;
                for my $alt (@ALL_AAS) {
                    unless ($ref eq $alt) {
                        print SUBS $ref.$pos.$alt."\n";
                    }
                }
            }

            close SUBS;
        }

        # and run sift on it

        $self->dbc->disconnect_when_inactive(1);

        my $cmd = "$sift_dir/bin/info_on_seqs $aln_file $subs_file $res_file";

        system($cmd) == 0 or die "Failed to run $cmd: $?";

        $self->dbc->disconnect_when_inactive(0);

        # parse and store the results 
        
        open (RESULTS, "<$res_file") or die "Failed to open $res_file: $!";

        # parse the results file

        my $pred_matrix = Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix->new(
            -analysis       => 'sift',
            -peptide_length => length($peptide),
        );

        while (<RESULTS>) {

            chomp;

            next if /WARNING/;
            next if /NOT SCORED/;

            my ($subst, $prediction, $score, $median_cons, $num_seqs, $blocks) = split;

            my ($ref_aa, $pos, $alt_aa) = $subst =~ /([A-Z])(\d+)([A-Z])/;

            next unless $ref_aa && $alt_aa && defined $pos;

            $pred_matrix->add_prediction(
                $pos,
                $alt_aa,
                $prediction, 
                $score,
            );
        }
        
        $self->save_predictions($pred_matrix);
    }

    # tar up the files

    system("tar --remove-files --exclude *.tgz -czvf $tarball * > /dev/null") == 0
        or die "Failed to create $output_dir/$tarball: $!";
}
