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
package Bio::EnsEMBL::Variation::Pipeline::ProteinFunction::RunSift;

use strict;

use File::Path qw(make_path remove_tree);
use Data::Dumper;

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix qw(@ALL_AAS);
use Bio::EnsEMBL::Variation::Utils::ComparaUtils qw(dump_alignment_for_sift);

use base ('Bio::EnsEMBL::Variation::Pipeline::ProteinFunction::BaseProteinFunction');

my $MEDIAN_CUTOFF = 2.75; # as per README

sub run {
  my $self = shift;

  my $translation_md5 = $self->required_param('translation_md5');
  my $sift_dir        = $self->required_param('sift_dir');
  my $working_dir     = $self->required_param('sift_working');
  my $ncbi_dir        = $self->required_param('ncbi_dir');
  my $blastdb         = $self->required_param('blastdb');

  my $dir = substr($translation_md5, 0, 2);
  my $output_dir = "$working_dir/$dir/$translation_md5";
  
  my $tarball = 'scratch.tgz';

  unless (-d $output_dir) {
    my $err;
    make_path($output_dir, {error => \$err});
    die "make_path failed: ".Dumper($err) if $err && @$err;
  }

  chdir $output_dir or die "Failed to chdir to $output_dir";

  my $fasta_file = "protein.fa";
  my $aln_file   = "protein.alignedfasta";
  my $res_file   = "protein.SIFTprediction";
  my $subs_file  = "subs.txt";

  if (-e "$output_dir/$tarball") {
    system("tar zxvf $tarball > /dev/null") == 0
      or die "Failed to untar $output_dir/$tarball: $!";
  }
  
  # set necessary environment variables for sift

  $ENV{NCBI}     = $ncbi_dir;
  $ENV{BLIMPS_DIR} = $sift_dir.'/blimps';
  $ENV{SIFT_DIR}   = $sift_dir;   
  $ENV{tmpdir}   = $output_dir;   
  
  # fetch our protein 

  my $peptide = $self->get_protein_sequence($translation_md5);

  my $alignment_ok = 1;

  unless (-e $aln_file) {
    
    # we need to get the multiple alignment
    
    if ($self->param('use_compara')) {

      my $stable_id = $self->get_stable_id_for_md5($translation_md5);
      
      eval {
        dump_alignment_for_sift($stable_id, $aln_file);
      };

      if ($@) {
        warn "Failed to get a compara alignment for $stable_id: $@";
        $alignment_ok = 0;
      }
    }
    else {

      # do the alignment ourselves
      
      # first create a fasta file for the protein sequence

      open (FASTA_FILE, ">$fasta_file");
      
      my $pep_copy = $peptide;
      $pep_copy =~ s/(.{80})/$1\n/g;
      chomp $pep_copy;
      print FASTA_FILE ">$translation_md5\n$pep_copy\n";

      close FASTA_FILE;

      # and run the alignment program

      $self->dbc->disconnect_when_inactive(1);

      my $cmd = "$sift_dir/bin/ensembl_seqs_chosen_via_median_info.csh $fasta_file $blastdb $MEDIAN_CUTOFF";
      my ($exit_code, $stderr, $flat_cmd) = $self->run_system_command($cmd);

      #die `env`."\n".$cmd;

      # my $exit_code = system($cmd);

      $self->dbc->disconnect_when_inactive(0);
      
      if ($exit_code == 0) {
        $alignment_ok = 1;
      }
      else {
        # the alignment failed for some reason, what to do?
        die "Alignment for $translation_md5 failed - cmd: $flat_cmd: $stderr";
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
    my ($exit_code, $stderr, $flat_cmd) = $self->run_system_command($cmd);

    die("Failed to run $flat_cmd: $stderr\n") unless $exit_code == 0;

    $self->dbc->disconnect_when_inactive(0);

    # parse and store the results 
    
    open (RESULTS, "<$res_file") or die "Failed to open $res_file: $!";

    # parse the results file

    my $pred_matrix = Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix->new(
      -analysis       => 'sift',
      -peptide_length   => length($peptide),
      -translation_md5  => $translation_md5,
    );

    my %evidence_stored;
    my $results_available = 0;

    while (<RESULTS>) {

      chomp;

      next if /WARNING/;
      next if /NOT SCORED/;

      my ($subst, $prediction, $score, $median_cons, $num_seqs, $blocks) = split;

      my ($ref_aa, $pos, $alt_aa) = $subst =~ /([A-Z])(\d+)([A-Z])/;

      next unless $ref_aa && $alt_aa && defined $pos;

      $results_available = 1;

      my $low_quality = 0;
      $low_quality = 1 if $median_cons > 3.25 || $num_seqs < 10;

      $pred_matrix->add_prediction(
        $pos,
        $alt_aa,
        $prediction, 
        $score,
        $low_quality
      );
      unless ($evidence_stored{$pos} ==1) {
        ## add attribs by position
        $pred_matrix->add_evidence( 'sequence_number',  $pos, $num_seqs );
        $pred_matrix->add_evidence( 'conservation_score', $pos, $median_cons  );
        $evidence_stored{$pos} = 1;
      }
    }
    if ($results_available == 1 ){  
      # avoid entering null matrices
      my $var_dba = $self->get_species_adaptor('variation');

      my $pfpma = $var_dba->get_ProteinFunctionPredictionMatrixAdaptor
        or die "Failed to get matrix adaptor";
    
      $pfpma->store($pred_matrix);
    }
  }

  # tar up the files
  my ($exit_code, $stderr, $flat_cmd) = $self->run_system_command(
    "tar --remove-files --exclude *.tgz -czvf $tarball * > /dev/null"
  );
  
  die "Failed to create $output_dir/$tarball: $stderr" unless $exit_code == 0;
}
