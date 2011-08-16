package Bio::EnsEMBL::Variation::Pipeline::ProteinFunction::RunPolyPhen;

use strict;

use Digest::MD5 qw(md5_hex);
use File::Path qw(make_path remove_tree);
use Data::Dumper;

use Bio::EnsEMBL::Variation::Utils::ProteinFunctionUtils qw(@ALL_AAS $AA_LOOKUP);
use Bio::EnsEMBL::Variation::Utils::ComparaUtils qw(dump_alignment_for_polyphen);

use base ('Bio::EnsEMBL::Variation::Pipeline::ProteinFunction::BaseProteinFunction');

my $DEBUG           = 0;
my $PERSIST         = 1;
my $MAX_PSIC_SEQS   = 8190;
my $MAX_PSIC_SEQLEN = 409650;

sub run {
    my $self = shift;

    my $translation_stable_id = $self->required_param('translation_stable_id'); 

    my $pph_dir = $self->required_param('pph_dir');
    
    my $md5 = md5_hex($translation_stable_id);
    my $dir = substr($md5, 0, 2);
    my $output_dir = "$pph_dir/working/$dir/$translation_stable_id";

    my @to_delete;

    my @pph_dirs    = qw{alignments blastfiles profiles structures lock source_files};
    my @pipe_dirs   = qw{features errors};
    
    my $tarball = 'scratch.tgz';
 
    if ( (-e $output_dir) && (!$PERSIST) ) {
        # if we're not being persistent delete the old directory first
        chdir $output_dir or die "Failed to chdir to $output_dir";
        my $err;
        remove_tree(@pph_dirs, @pipe_dirs, {error => \$err});
        die "remove_tree failed: ".Dumper($err) if $err && @$err;
    }

    unless (-d $output_dir) {
        my $err;
        make_path($output_dir, {error => \$err});
        die "make_path failed: ".Dumper($err) if $err && @$err;
    }

    chdir $output_dir or die "Failed to chdir to $output_dir";
   
    if (-e "$output_dir/$tarball") {
        system("tar zxvf $tarball > /dev/null") == 0
            or die "Failed to untar $output_dir/$tarball: $!";
    }
    else {
        my $err;
        make_path(@pph_dirs, @pipe_dirs, {error => \$err});
        die "make_path failed: ".Dumper($err) if $err && @$err;
    }

    my $subs_file       = "${output_dir}/source_files/${translation_stable_id}_subs.txt";
    my $protein_file    = "${output_dir}/source_files/${translation_stable_id}_protein.fa";
    my $aln_file        = "${output_dir}/alignments/${translation_stable_id}.aln";
    my $output_file     = "${output_dir}/features/${translation_stable_id}.features";
    my $error_file      = "${output_dir}/errors/${translation_stable_id}.polyphen_stderr";

    my $tfa = $self->get_transcript_file_adaptor;

    # dump the protein sequence,

    open (PROTEIN, ">$protein_file") or die "Failed to open file for protein $protein_file: $!";
    
    print PROTEIN $tfa->get_translation_fasta($translation_stable_id);

    # and the substitutions.

    open (SUBS, ">$subs_file") or die "Failed to open file for SUBS $subs_file: $!";

    #push @to_delete, $subs_file, $protein_file;

    my $peptide = $tfa->get_translation_seq($translation_stable_id);
        
    die "$translation_stable_id is length 0?" unless length($peptide) > 0;

    my @aas = split //, $peptide;

    my $idx = 0;

    for my $ref (@aas) {
        $idx++;
         
        # ignore any non standard amino acids, e.g. X

        next unless defined $AA_LOOKUP->{$ref};        
        
        for my $alt (@ALL_AAS) {
            unless ($ref eq $alt) {
                print SUBS join ("\t",
                    $translation_stable_id.'_'.$idx.'_'.$alt,
                    $translation_stable_id,
                    $idx,
                    $ref,
                    $alt
                ), "\n";
            }
        }
    }

    close SUBS;
    close PROTEIN;
    
    if ($self->param('use_compara')) {
        
        # if we're using compara alignments then dump the alignment
        # to the alignment file, PolyPhen will check if it exists
        # and use it in place of its own alignment if so
        eval {
            dump_alignment_for_polyphen($translation_stable_id, $aln_file);
        };

        if ($@) {
            die "Failed to fetch a compara alignment for $translation_stable_id: $@";
        }
    }

    # now run polyphen itself, disconnecting for the duration

    $self->dbc->disconnect_when_inactive(1);

    my $cmd = "$pph_dir/bin/run_pph.pl -d $output_dir -s $protein_file $subs_file 1> $output_file 2> $error_file";

    system($cmd) == 0 or die "Failed to run $cmd: $?";
    
    $self->dbc->disconnect_when_inactive(0);
    
    my $scratch_dirs = join ' ', @pph_dirs;

    system("tar czvf $tarball $scratch_dirs > /dev/null") == 0 
        or die "tar command failed: $?";
    
    my $err;
    remove_tree(@pph_dirs, {error => \$err});
    die "remove_tree failed: ".Dumper($err) if $err && @$err;

    my $exit_code = system("gzip -f $output_file");
    
    if ($exit_code == 0) {
        $output_file .= '.gz';
    }
    else {
         warn "Failed to gzip $output_file: $?"; 
    }

    $self->param('feature_file', $output_file);

    if (-s $error_file) {
        warn "run_pph.pl STDERR output in $error_file\n";
    }
    else {
        push @to_delete, $error_file;
    }

    # delete unnecesary files

    unlink @to_delete;
}

sub write_output {
    my $self = shift;
    
    if (my $feature_file = $self->param('feature_file')) {
        $self->dataflow_output_id( [{
            translation_stable_id   => $self->param('translation_stable_id'),
            feature_file            => $feature_file,
        }], 2);
    }
}

1;

