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

use base ('Bio::EnsEMBL::Variation::Pipeline::ProteinFunction::BaseProteinFunction');

my $CADD_CUTOFF = 0.0;
my $REVEL_CUTOFF = 0.5;
sub run {
  my $self = shift;
  my $working_dir = $self->required_param('dbnsfp_working');
  my $dbnsfp_file = $self->required_param('dbnsfp_file');
  unless (-d $working_dir) {
    my $err;
    make_path($working_dir, {error => \$err});
    die "make_path failed: ".Dumper($err) if $err && @$err;
  }

  my $codonTable = Bio::Tools::CodonTable->new();

  my $predictions = {
    dbnsfp_meta_svm => {
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

  my $vdba = $self->get_species_adaptor('variation');
  my $cdba = $self->get_species_adaptor('core');

  my $slice_adaptor = $cdba->get_SliceAdaptor or die "Failed to get slice adaptor";
  my $translation_adaptor = $cdba->get_TranslationAdaptor or die "Failed to get translation adaptor";
  my $pfpma = $vdba->get_ProteinFunctionPredictionMatrixAdaptor or die "Failed to get matrix adaptor";

  my $translation_md5 = $self->required_param('translation_md5');
  my $translation_stable_id = $self->get_stable_id_for_md5($translation_md5);
  my $translation = $translation_adaptor->fetch_by_stable_id($translation_stable_id);

  my $pred_matrices = {};
  my $results_available = {};
  foreach my $analysis (qw/dbnsfp_cadd dbnsfp_revel dbnsfp_meta_svm dbnsfp_mutation_assessor/) {
    my $pred_matrix = Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix->new(
      -analysis       => $analysis,
      -peptide_length   => $translation->length,
      -translation_md5  => $translation_md5,
    );
    $pred_matrices->{$analysis} = $pred_matrix;
  }

  my $transcript = $translation->transcript;
  my $chrom = $transcript->seq_region_name;
  my $start = $transcript->seq_region_start;
  my $end = $transcript->seq_region_end;
  my $strand = $transcript->seq_region_strand;
  my $slice = $slice_adaptor->fetch_by_region('chromosome', $chrom,  $start, $end);
  my $transcript_mapper = $transcript->get_TranscriptMapper();

  my $transcript_stable_id =  $transcript->stable_id;
  my $translation_seq = $translation->seq;
  my @amino_acids = ();
  foreach my $i (1 .. $translation->length) {
    my @pep_coordinates = $transcript_mapper->pep2genomic($i, $i);
    my $triplet = '';
    my @coords = ();
    foreach my $coord (@pep_coordinates) {
      my $coord_start = $coord->start;
      my $coord_end = $coord->end;
      next if ($coord_start <= 0);
      my $new_start = $coord_start - $start + 1;
      my $new_end   = $coord_end   - $start + 1;
      my $subseq = $slice->subseq($new_start, $new_end, $strand);
      $triplet .= $subseq;
      push @coords, [$coord_start, $coord_end];
    }
    my $aa = $codonTable->translate($triplet);
    if (!$aa) {
      push @amino_acids, 'X';
      next;
    }
    push @amino_acids, $aa;
    
    my $new_triplets = mutate($triplet);
    foreach my $coord (@coords) {
      my $triplet_start = $coord->[0];
      my $triplet_end = $coord->[1];
      my $iter = $obj->query("$chrom:$triplet_start-$triplet_end");
      while (my $line = $iter->next) {
        $line =~ s/\r$//g;
        my @split = split /\t/, $line;
        my %data = map {$headers->[$_] => $split[$_]} (0..(scalar @{$headers} - 1));
        my $chr = $data{'#chr'};
        my $cadd_raw = $data{'CADD_raw'};
        my $cadd_rank = $data{'CADD_raw_rankscore'};
        my $cadd_phred = $data{'CADD_phred'};
        my $revel_raw = $data{'REVEL_score'};
        my $revel_rank = $data{'REVEL_rankscore'};
        my $metaSVM_score = $data{'MetaSVM_score'};
        my $metaSVM_pred = $data{'MetaSVM_pred'};
        my $mutation_assessor_score = $data{'MutationAssessor_score'};
        my $mutation_assessor_pred = $data{'MutationAssessor_pred'};
        my $pos = $data{'pos(1-based)'};
        my $nucleotide_position = $pos - $triplet_start;
        my $ref = $data{'ref'};
        my $refcodon = $data{'refcodon'};
        my $alt = $data{'alt'};
        next if ($alt eq $ref);
        my $aaalt = $data{'aaalt'};
        my $aaref = $data{'aaref'};
        my $mutated_triplet =  $new_triplets->{$triplet}->{$nucleotide_position}->{$alt};
        my $mutated_aa = $codonTable->translate($mutated_triplet);
        if ($cadd_raw ne '.') {
          $results_available->{'dbnsfp_cadd'} = 1; 
          my $prediction = ($cadd_raw >= $CADD_CUTOFF) ? 'simulated' : 'observed';
          my $low_quality = 0;
          $pred_matrices->{'dbnsfp_cadd'}->add_prediction(
            $i,
            $mutated_aa,
            $prediction,
            $cadd_raw,
            $low_quality,
          );
        }
        if ($revel_raw ne '.') {
          $results_available->{'dbnsfp_revel'} = 1; 
          my $prediction = ($revel_raw >= $REVEL_CUTOFF) ? 'likely_disease_causing' : 'likely_not_disease_causing';
          my $low_quality = 0;
          $pred_matrices->{'dbnsfp_revel'}->add_prediction(
            $i,
            $mutated_aa,
            $prediction,
            $revel_raw,
            $low_quality,
          );
        }
        if ($metaSVM_score ne '.') {
          $results_available->{'dbnsfp_meta_svm'} = 1; 
          my $prediction = $predictions->{dbnsfp_meta_svm}->{$metaSVM_pred};
          my $low_quality = 0;
          $pred_matrices->{'dbnsfp_meta_svm'}->add_prediction(
            $i,
            $mutated_aa,
            $prediction,
            $metaSVM_score,
            $low_quality,
          );
        }
        if ($mutation_assessor_score ne '.') {
          $results_available->{'dbnsfp_mutation_assessor'} = 1; 
          my $prediction = $predictions->{dbnsfp_mutation_assessor}->{$mutation_assessor_pred};
          my $low_quality = 0;
          $pred_matrices->{'dbnsfp_mutation_assessor'}->add_prediction(
            $i,
            $mutated_aa,
            $prediction,
            $mutation_assessor_score,
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
    }
  }
}

sub mutate {
  my $triplet = shift;
  my @nucleotides = split('', $triplet);
  my $new_triplets;
  foreach my $i (0 .. $#nucleotides) {
    my $mutations = get_mutations();
    $new_triplets = get_mutated_triplets($triplet, $mutations, $i, $new_triplets);
  }
  return $new_triplets;
}

sub get_mutated_triplets {
  my $triplet = shift;
  my $mutations = shift;
  my $position = shift;
  my $new_triplets = shift;
  foreach my $mutation (@$mutations) {
    my $update_triplet = $triplet;
    substr($update_triplet, $position, 1, $mutation);
    $new_triplets->{$triplet}->{$position}->{$mutation} = $update_triplet;
  }
  return $new_triplets;
}

sub get_mutations {
  return ['A', 'G', 'C', 'T'];
}

1;
