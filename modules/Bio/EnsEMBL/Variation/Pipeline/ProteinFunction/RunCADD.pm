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
package Bio::EnsEMBL::Variation::Pipeline::ProteinFunction::RunCADD;

use strict;
use Bio::DB::HTS::Tabix;
use Bio::Tools::CodonTable;
use File::Path qw(make_path remove_tree);
use FileHandle;
use Data::Dumper;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix qw(@ALL_AAS);

use base ('Bio::EnsEMBL::Variation::Pipeline::ProteinFunction::BaseProteinFunction');

my $CADD_CUTOFF = 30;

sub run {
  my $self = shift;
  $self->dbc and $self->dbc->disconnect_if_idle();
  $self->dbc->disconnect_when_inactive(1);

  my $working_dir = $self->required_param('cadd_working');
  my $cadd_file = $self->required_param('cadd_file');
  unless (-d $working_dir) {
    my $err;
    make_path($working_dir, {error => \$err});
    die "make_path failed: ".Dumper($err) if $err && @$err;
  }

  my $codonTable = Bio::Tools::CodonTable->new();

  my $obj = Bio::DB::HTS::Tabix->new(filename => $cadd_file);

  my $headers;
  open HEAD, "tabix -fh $cadd_file 1:1-1 2>&1 | ";
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

  my $results_available = 0;
  my $pred_matrix = Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix->new(
    -analysis        => 'cadd',
    -peptide_length  => $translation->length,
    -translation_md5 => $translation_md5,
  );

  my @amino_acids = ();
  my @all_triplets = @{$self->get_triplets($translation_stable_id)};

  my $debug_data = {};

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
        my $chr = $data{'#Chr'};
        my $pos = $data{'Pos'};
        my $ref = $data{'Ref'};
        my $alt = $data{'Alt'};
        my $cadd_phred = $data{'PHRED'};
        next if ($alt eq $ref);
        my $nucleotide_position = ($reverse) ? $triplet_end - $pos : $pos - $triplet_start;
        my $mutated_triplet = $new_triplets->{$triplet_seq}->{$nucleotide_position}->{$alt};
        my $mutated_aa = $codonTable->translate($mutated_triplet);

        if ($cadd_phred ne '.') {
          $results_available = 1;
          my $prediction = ($cadd_phred >= $CADD_CUTOFF) ? 'likely deleterious' : 'likely benign';
          my $low_quality = 0;
          if ($self->param('debug_mode')) {
            $debug_data->{$i}->{$mutated_aa}->{$prediction} = $cadd_phred;
          }
          $pred_matrix->add_prediction(
            $i,
            $mutated_aa,
            $prediction,
            $cadd_phred,
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
  if ($results_available) {
    $pfpma->store($pred_matrix);
  }
  if ($self->param('debug_mode')) {
    my $matrix = $pfpma->fetch_by_analysis_translation_md5('cadd', $translation_md5);
    my $fh = FileHandle->new("$working_dir/debug_$translation_stable_id", 'w');
    foreach my $i (keys %$debug_data) {
      foreach my $aa (keys %{$debug_data->{$i}}) {
        next if ($aa eq '*');
        foreach my $prediction (keys %{$debug_data->{$i}->{$aa}}) {
          my ($new_pred, $new_score) = $matrix->get_prediction($i, $aa); 
          print $fh join(' ', $i, $aa, $prediction, $debug_data->{$i}->{$aa}->{$prediction}, $new_pred, $new_score), "\n";
        }
      }
    }
    $fh->close;
  }

}

1;
