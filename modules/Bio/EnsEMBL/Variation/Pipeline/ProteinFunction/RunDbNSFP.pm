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
package Bio::EnsEMBL::Variation::Pipeline::ProteinFunction::RunDbNSFP;

use strict;
use Bio::DB::HTS::Tabix;
use Bio::Tools::CodonTable;
use File::Path qw(make_path remove_tree);
use Data::Dumper;
use FileHandle;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix qw(@ALL_AAS);
use Data::Dumper;

use base ('Bio::EnsEMBL::Variation::Pipeline::ProteinFunction::BaseProteinFunction');

my $REVEL_CUTOFF = 0.5;

sub run {
  my $self = shift;
  $self->dbc and $self->dbc->disconnect_if_idle();
  $self->dbc->disconnect_when_inactive(1);

  my $working_dir = $self->required_param('dbnsfp_working');
  my $dbnsfp_file = $self->required_param('dbnsfp_file');
  unless (-d $working_dir) {
    my $err;
    make_path($working_dir, {error => \$err});
    die "make_path failed: ".Dumper($err) if $err && @$err;
  }

  my $codonTable = Bio::Tools::CodonTable->new();

  my $predictions = {
    dbnsfp_meta_lr => {
      T => 'tolerated',
      D => 'damaging',
    },
    dbnsfp_mutation_assessor => {
      H => 'high',
      M => 'medium',
      L => 'low',
      N => 'neutral',
    }
  };

  my $obj = Bio::DB::HTS::Tabix->new(filename => $dbnsfp_file);

  my $headers;
  open HEAD, "tabix -fh $dbnsfp_file 1:1-1 2>&1 | ";
  while(<HEAD>) {
    next unless /^\#/;
    chomp;
    $headers = [split];
  }
  close HEAD;

  my $translation_md5 = $self->required_param('translation_md5');
  my $translation_stable_id = $self->get_stable_id_for_md5($translation_md5);

  if (!$translation_stable_id) {
    $self->warning("No stable id for $translation_md5");
    return;
  }

  my $translation = $self->get_translation($translation_stable_id);
  my $translation_seq = $translation->seq;
  my $transcript = $translation->transcript;
  my $reverse = $transcript->strand < 0;
  my $transcript_stable_id = $transcript->stable_id;

  my $vdba = $self->get_species_adaptor('variation');
  my $pfpma = $vdba->get_ProteinFunctionPredictionMatrixAdaptor or die "Failed to get matrix adaptor";

  my $pred_matrices = {};
  my $results_available = {};
  foreach my $analysis (qw/dbnsfp_revel dbnsfp_meta_lr dbnsfp_mutation_assessor/) {
    my $pred_matrix = Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix->new(
      -analysis       => $analysis,
      -peptide_length   => $translation->length,
      -translation_md5  => $translation_md5,
    );
    $pred_matrices->{$analysis} = $pred_matrix;
  }

  my $debug_data = {};

  my @amino_acids = ();
  my @all_triplets = @{$self->get_triplets($translation_stable_id)};
  foreach my $entry (@all_triplets) {
    my $aa = $entry->{aa};
    push @amino_acids, $aa;
    next if $aa eq 'X';
    my @coords = @{$entry->{coords}};
    my $chrom = $entry->{chrom};
    my $triplet_seq = $entry->{triplet_seq};
    my $i = $entry->{aa_position};
    my $new_triplets = $entry->{new_triplets};
    foreach my $coord (@coords) {
      my $triplet_start = $coord->[0];
      my $triplet_end = $coord->[1];
      my $iter = $obj->query("$chrom:$triplet_start-$triplet_end");
      while (my $line = $iter->next) {
        $line =~ s/\r$//g;
        my @split = split /\t/, $line;
        my %data = map {$headers->[$_] => $split[$_]} (0..(scalar @{$headers} - 1));
        my $chr = $data{'#chr'};
        my $revel_raw = $data{'REVEL_score'};
        my $revel_rank = $data{'REVEL_rankscore'};
        my $metaLR_score = $data{'MetaLR_score'};
        my $metaLR_pred = $data{'MetaLR_pred'};
        my $mutation_assessor_rankscore = $data{'MutationAssessor_score_rankscore'};
        my $mutation_assessor_pred = $data{'MutationAssessor_pred'};
        my $pos = $data{'pos(1-based)'};
        my $nucleotide_position = ($reverse) ? $triplet_end - $pos : $pos - $triplet_start;
        my $ref = $data{'ref'};
        my $refcodon = $data{'refcodon'};
        my $alt = $data{'alt'};
        next if ($alt eq $ref);
        my $aaalt = $data{'aaalt'};
        my $aaref = $data{'aaref'};
        my $mutated_triplet =  $new_triplets->{$triplet_seq}->{$nucleotide_position}->{$alt};
        my $mutated_aa = $codonTable->translate($mutated_triplet);
        next if ($aaalt ne $mutated_aa);
        if ($revel_raw ne '.') {
          $results_available->{'dbnsfp_revel'} = 1; 
          my $prediction = ($revel_raw >= $REVEL_CUTOFF) ? 'likely disease causing' : 'likely benign';
          my $low_quality = 0;
          $debug_data->{dbnsfp_revel}->{$i}->{$mutated_aa}->{$prediction} = $revel_raw if ($self->param('debug_mode'));
          $pred_matrices->{'dbnsfp_revel'}->add_prediction(
            $i,
            $mutated_aa,
            $prediction,
            $revel_raw,
            $low_quality,
          );
        }
        if ($metaLR_score ne '.') {
          $results_available->{'dbnsfp_meta_lr'} = 1; 
          my $prediction = $predictions->{dbnsfp_meta_lr}->{$metaLR_pred};
          my $low_quality = 0;
          $debug_data->{dbnsfp_meta_lr}->{$i}->{$mutated_aa}->{$prediction} = $metaLR_score if ($self->param('debug_mode'));
          $pred_matrices->{'dbnsfp_meta_lr'}->add_prediction(
            $i,
            $mutated_aa,
            $prediction,
            $metaLR_score,
            $low_quality,
          );
        }
        if ($mutation_assessor_rankscore ne '.') {
          $results_available->{'dbnsfp_mutation_assessor'} = 1; 
          my $prediction = $predictions->{dbnsfp_mutation_assessor}->{$mutation_assessor_pred};
          my $low_quality = 0;
          $debug_data->{dbnsfp_mutation_assessor}->{$i}->{$mutated_aa}->{$prediction} = $mutation_assessor_rankscore if ($self->param('debug_mode'));
          $pred_matrices->{'dbnsfp_mutation_assessor'}->add_prediction(
            $i,
            $mutated_aa,
            $prediction,
            $mutation_assessor_rankscore,
            $low_quality,
          );
        }
      }
    }
  }

  if ($translation_seq ne join('', @amino_acids)) {
    my $fh = FileHandle->new("$working_dir/$translation_stable_id", 'w');
    print $fh "$transcript_stable_id\n$translation_seq\n";
    print $fh join('', @amino_acids), "\n";
    $fh->close;
  }
  foreach my $analysis (keys %$pred_matrices) {
    my $pred_matrix = $pred_matrices->{$analysis};
    if ($results_available->{$analysis}) {
      $pfpma->store($pred_matrix);
      if ($self->param('debug_mode')) {
        my $fh = FileHandle->new("$working_dir/$analysis\_$translation_stable_id", 'w');
        my $matrix = $pfpma->fetch_by_analysis_translation_md5($analysis, $translation_md5); 
        foreach my $i (keys %{$debug_data->{$analysis}}) {
          foreach my $aa (keys %{$debug_data->{$analysis}->{$i}}) {
            next if ($aa eq '*');
            foreach my $prediction (keys %{$debug_data->{$analysis}->{$i}->{$aa}}) {
              my ($new_pred, $new_score) = $matrix->get_prediction($i, $aa);
              print $fh join(' ', $analysis, $i, $aa, $prediction, $debug_data->{$analysis}->{$i}->{$aa}->{$prediction}, $new_pred, $new_score), "\n";
            }
          }
        }
        $fh->close;
      }
    }
  }
}

1;
