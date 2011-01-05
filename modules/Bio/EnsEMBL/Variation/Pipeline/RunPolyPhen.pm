=head1 LICENSE

 Copyright (c) 1999-2011 The European Bioinformatics Institute and
 Genome Research Limited.  All rights reserved.

 This software is distributed under a modified Apache license.
 For license details, please see

   http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <dev@ensembl.org>.

 Questions may also be sent to the Ensembl help desk at
 <helpdesk@ensembl.org>.

=cut

package Bio::EnsEMBL::Hive::RunnableDB::nsSNP::RunPolyPhen;

use strict;

use Digest::MD5 qw(md5_hex);
use File::Path qw(make_path remove_tree);
use Data::Dumper;

#use Bio::EnsEMBL::Registry;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;

use base ('Bio::EnsEMBL::Hive::RunnableDB::nsSNP::BaseRunnable');

my $DEBUG   = 0;

sub _run {
    my $self = shift;

    $self->dbc->disconnect_when_inactive(1);
    
    my $gene_id = $self->param('gene_stable_id') 
        or die "gene_stable_id is a required parameter";

    my $pph_dir = $self->param('pph_dir')
        or die "pph_dir is a required parameter";

#    my $reg = 'Bio::EnsEMBL::Registry';
#    $reg->load_all;
#    my $ga = $reg->get_adaptor('human', 'core', 'gene') 
#        or die "Failed to get gene adaptor";
#    my $tva = $reg->get_adaptor('human', 'variation', 'transcriptvariation')
#        or die "Failed to get transcript variation adaptor";
#    $reg->set_disconnect_when_inactive;
 
    my $core_dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new( 
        '-species' => 'Homo_sapiens',
        '-group'   => 'core',
        '-port'    => 3306,
        '-host'    => 'ens-staging',
        '-user'    => 'ensro',
        '-pass'    => '',
        '-dbname'  => 'homo_sapiens_core_59_37d',
    );

    my $var_dba = Bio::EnsEMBL::Variation::DBSQL::DBAdaptor->new(
        '-species' => 'Homo_sapiens',
        '-group'   => 'variation',
        '-port'    => 3306,
        '-host'    => 'ens-staging',
        '-user'    => 'ensro',
        '-pass'    => '',
        '-dbname'  => 'homo_sapiens_variation_59_37d', 
    );
   
    my $ga = $core_dba->get_GeneAdaptor 
        or die "Failed to get gene adaptor";

    my $tva = $var_dba->get_TranscriptVariationAdaptor
        or die "Failed to get transcript variation adaptor";
   
    my $gene = $ga->fetch_by_stable_id($gene_id) 
        or die "Failed to fetch gene for stable id: $gene_id";

    my %proteins;

    # create a directory tree to control number of files in a directory

    my $md5 = md5_hex($gene_id);
    my $dir = substr($md5, 0, 2);
    #$dir .= '/'.substr($md5, 2, 2);
    my $output_dir = "$pph_dir/working/$dir/$gene_id";

    my @pph_dirs    = qw{alignments blastfiles profiles structures lock};
    my @pipe_dirs   = qw{features errors};
    
    my $tarball = 'scratch.tgz';
  
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

    my $snp_file        = "/tmp/${gene_id}_snps.txt";
    my $proteins_file   = "/tmp/${gene_id}_proteins.fa";
    my $output_file     = "${output_dir}/features/${gene_id}.features";
    my $error_file      = "${output_dir}/errors/${gene_id}.polyphen_stderr";

    open (SNPS, ">$snp_file") or die "Failed to open file for SNPs: $!";
    open (PROTEINS, ">$proteins_file") or die "Failed to open file for proteins: $!";

    $self->delete_after($snp_file, $proteins_file);

    my $found_snp = 0;

    for my $tran (@{ $gene->get_all_Transcripts }) {
        for my $tv (@{ $tva->fetch_all_by_Transcripts([$tran]) }) {
            if ($tv->consequence_type->[0] eq 'NON_SYNONYMOUS_CODING') {
                my $tl = $tran->translation;
                my ($aa1, $aa2) = split /\//, $tv->pep_allele_string;

                print SNPS join ("\t",
                    $tv->variation_feature->variation->name,
                    $tl->stable_id,
                    $tv->translation_start,
                    $aa1,
                    $aa2
                ), "\n";

                unless ($proteins{$tl->stable_id}++) {
                    my $seq = $tl->seq;
                    $seq =~ s/(.{80})/$1\n/g;
                    print PROTEINS '>'.$tl->stable_id."\n$seq\n";
                }

                $found_snp = 1;
            }
        }
    }

    close SNPS;
    close PROTEINS;
   

    if ($found_snp) {

        my $cmd = "$pph_dir/bin/run_pph.pl -d $output_dir -s $proteins_file $snp_file 1> $output_file 2> $error_file";

        system($cmd) == 0 or die "Failed to run $cmd: $?";
        
        my $scratch_dirs = join ' ', @pph_dirs;

        system("tar czvf scratch.tgz $scratch_dirs > /dev/null") == 0 
            or die "tar command failed: $?";
        
        my $err;
        remove_tree(@pph_dirs, {error => \$err});
        die "remove_tree failed: ".Dumper($err) if $err && @$err;

        $self->param('feature_file', $output_file);

        if (-s $error_file) {
            warn "run_pph.pl STDERR output in $error_file\n";
        }
        else {
            $self->delete_after($error_file);
        }
    }
    else {
        my $err;
        remove_tree(@pph_dirs, @pipe_dirs, {error => \$err});
        die "remove_tree failed: ".Dumper($err) if $err && @$err;
    }
}

sub _write_output {
    my $self = shift;
    
    if (my $feature_file = $self->param('feature_file')) {
        $self->dataflow_output_id( {
            gene_stable_id  => $self->param('gene_stable_id'),
            pph_dir         => $self->param('pph_dir'),
            feature_file    => $feature_file
        }, 3);
    }
}

1;
