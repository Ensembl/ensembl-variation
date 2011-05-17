package Bio::EnsEMBL::Variation::Pipeline::ProteinFunction::RunPolyPhen;

use strict;

use Digest::MD5 qw(md5_hex);
use File::Path qw(make_path remove_tree);
use Data::Dumper;

use base ('Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess');

my $DEBUG   = 0;

my $PERSIST = 0;

sub run {
    my $self = shift;

    my $transcript_stable_id = $self->param('transcript_stable_id'); 

    my $pph_dir = $self->param('pph_dir');
 
    my @all_aas = qw(A C D E F G H I K L M N P Q R S T V W Y);
    
    my $aa_muts = {
        A => [qw(D E G P S T V)],
        C => [qw(F G R S W Y)],
        D => [qw(A E G H N V Y)],
        E => [qw(A D G K Q V)],
        F => [qw(C I L S V Y)],
        G => [qw(A C D E R S V W)],
        H => [qw(D L N P Q R Y)],
        I => [qw(F K L M N R S T V)],
        K => [qw(E I M N Q R T)],
        L => [qw(F H I M P Q R S V W)],
        M => [qw(I K L R T V)],
        N => [qw(D H I K S T Y)],
        P => [qw(A H L Q R S T)],
        Q => [qw(E H K L P R)],
        R => [qw(C G H I K L M P Q S T W)],
        S => [qw(A C F G I L N P R T W Y)],
        T => [qw(A I K M N P R S)],
        V => [qw(A D E F G I L M)],
        W => [qw(C G L R S)],
        Y => [qw(C D F H N S)],
    };

    # create a directory tree to control number of files in a directory

    my $md5 = md5_hex($transcript_stable_id);
    my $dir = substr($md5, 0, 2);
    #$dir .= '/'.substr($md5, 2, 2);
    my $output_dir = "$pph_dir/working/$dir/$transcript_stable_id";

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

    my $snp_file        = "${output_dir}/source_files/${transcript_stable_id}_snps.txt";
    my $protein_file    = "${output_dir}/source_files/${transcript_stable_id}_protein.fa";
    my $output_file     = "${output_dir}/features/${transcript_stable_id}.features";
    my $error_file      = "${output_dir}/errors/${transcript_stable_id}.polyphen_stderr";

    open (SNPS, ">$snp_file") or die "Failed to open file for SNPs $snp_file: $!";
    open (PROTEIN, ">$protein_file") or die "Failed to open file for protein $protein_file: $!";
  
    #push @to_delete, $snp_file, $protein_file;

    my $peptide = $self->get_transcript_file_adaptor->get_protein_seq($transcript_stable_id);
        
    die "$transcript_stable_id has no peptide?" unless length($peptide) > 0;

    my @aas = split //, $peptide;

    my $tran_ver = $transcript_stable_id;

    $peptide =~ s/(.{80})/$1\n/g;
    print PROTEIN ">$tran_ver\n$peptide\n";

    my $idx = 0;

    for my $ref (@aas) {
        $idx++;
        #my $alt_aas = $aa_muts->{$aa};
        #next unless $alt_aas;
        next unless $aa_muts->{$ref}; # ignore any non standard amino acids, e.g. X
        for my $alt (@all_aas) {
            unless ($ref eq $alt) {
                print SNPS join ("\t",
                    $tran_ver.'_'.$idx.'_'.$alt,
                    $tran_ver,
                    $idx,
                    $ref,
                    $alt
                ), "\n";
            }
        }
    }

    close SNPS;
    close PROTEIN;

    $self->dbc->disconnect_when_inactive(1);

    my $cmd = "$pph_dir/bin/run_pph.pl -d $output_dir -s $protein_file $snp_file 1> $output_file 2> $error_file";

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
            transcript_stable_id    => $self->param('transcript_stable_id'),
            feature_file            => $feature_file,
        }], 2);
    }
}

1;

