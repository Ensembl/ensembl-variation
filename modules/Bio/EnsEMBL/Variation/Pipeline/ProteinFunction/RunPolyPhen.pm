=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut


=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

 Questions may also be sent to the Ensembl help desk at
 <http://www.ensembl.org/Help/Contact>.

=cut
package Bio::EnsEMBL::Variation::Pipeline::ProteinFunction::RunPolyPhen;

use strict;

use Digest::MD5 qw(md5_hex);
use File::Path qw(make_path remove_tree);
use Data::Dumper;

use Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix qw(@ALL_AAS $AA_LOOKUP);
use Bio::EnsEMBL::Variation::Utils::ComparaUtils qw(dump_alignment_for_polyphen);

use base ('Bio::EnsEMBL::Variation::Pipeline::ProteinFunction::BaseProteinFunction');

my $PERSIST         = 1;
my $MAX_PSIC_SEQS   = 8190;
my $MAX_PSIC_SEQLEN = 409650;

sub run {
    my $self = shift;

    my $translation_md5     = $self->required_param('translation_md5'); 

    my $pph_dir     = $self->required_param('pph_dir');
    my $working_dir = $self->required_param('pph_working');
    
    my $dir = substr($translation_md5, 0, 2);
    my $output_dir = "$working_dir/$dir/$translation_md5";

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

    my $subs_file       = "${output_dir}/source_files/subs.txt";
    my $protein_file    = "${output_dir}/source_files/protein.fa";
    my $aln_file        = "${output_dir}/alignments/${translation_md5}.aln";
    my $output_file     = "${output_dir}/features/features.txt";
    my $error_file      = "${output_dir}/errors/polyphen.err";

    # dump the protein sequence,

    my $peptide = $self->get_protein_sequence($translation_md5);

    open (PROTEIN, ">$protein_file") or die "Failed to open file for protein $protein_file: $!";
   
    my $pep_copy = $peptide;
    $pep_copy =~ s/(.{80})/$1\n/g;
    chomp $pep_copy;
    print PROTEIN ">$translation_md5\n$pep_copy\n";
   
    close PROTEIN;

    # and the substitutions.

    open (SUBS, ">$subs_file") or die "Failed to open file for SUBS $subs_file: $!";

    #push @to_delete, $subs_file, $protein_file;

    my @aas = split //, $peptide;

    my $idx = 0;

    for my $ref (@aas) {
        $idx++;
         
        # ignore any non standard amino acids, e.g. X

        next unless defined $AA_LOOKUP->{$ref};        
        
        for my $alt (@ALL_AAS) {
            unless ($ref eq $alt) {
                print SUBS join ("\t",
                    $translation_md5,
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
        
        my $stable_id = $self->get_stable_id_for_md5($translation_md5);
        
        # if we're using compara alignments then dump the alignment
        # to the alignment file, PolyPhen will check if it exists
        # and use it in place of its own alignment if so

        eval {
            dump_alignment_for_polyphen($stable_id, $aln_file);
        };

        if ($@) {
            die "Failed to fetch a compara alignment for $stable_id: $@";
        }
    }

    # now run polyphen itself, disconnecting for the duration

    $self->dbc->disconnect_when_inactive(1);

    # use -A option to disable polyphen's own LSF support (which conflicts with the hive)
    my $cmd = "$pph_dir/bin/run_pph.pl -A -d $output_dir -s $protein_file $subs_file 1> $output_file 2> $error_file";

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
            translation_md5 => $self->param('translation_md5'),
            feature_file    => $feature_file,
        }], 2);
    }
}

1;

