=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2021] EMBL-European Bioinformatics Institute

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
package Bio::EnsEMBL::Variation::Pipeline::ProteinFunction::BaseProteinFunction;

use strict;
use warnings;

use Bio::SimpleAlign;
use Bio::LocatableSeq;
use Bio::Tools::CodonTable;
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);

use base ('Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess');

sub get_protein_sequence {
    my ($self, $md5) = @_;

    my $tfa = $self->get_transcript_file_adaptor;

    my $fasta = $tfa->get_translation_fasta($md5);

    unless (length($fasta) > 34) {
        # sleep in case it's some race condition and then try again
        sleep(rand(5));
        $fasta = $tfa->get_translation_fasta($md5);

        unless (length($fasta) > 34) {
            die "$md5 looks weirdly short!";
        }
    }

    # strip out fasta header etc.

    $fasta =~ s/>.*\n//m;
    $fasta =~ s/\s//mg;
    
    die "No peptide for $md5?" unless length($fasta) > 0;

    return $fasta;
}

sub get_stable_id_for_md5 {
    my ($self, $md5) = @_;

    my $var_dba = $self->get_species_adaptor('variation');
    
    my $get_stable_id_sth = $var_dba->dbc->prepare(qq{
        SELECT  stable_id
        FROM    translation_mapping
        WHERE   md5 = ?
    });

    $get_stable_id_sth->execute($md5);

    my ($stable_id) = $get_stable_id_sth->fetchrow_array;

    return $stable_id;
}

sub get_translation {
  my $self = shift;
  my $translation_stable_id = shift;
  my $core_type = ($translation_stable_id =~ /^NP|XP/) ? 'otherfeatures' : 'core';
  my $cdba = $self->get_species_adaptor($core_type);
  my $translation_adaptor = $cdba->get_TranslationAdaptor or die "Failed to get translation adaptor";
  my $translation = $translation_adaptor->fetch_by_stable_id($translation_stable_id);
  return $translation;
}

sub mutate {
  my $self = shift;
  my $triplet = shift;
  my $reverse = shift;
  my @nucleotides = split('', $triplet);
  my $new_triplets;
  foreach my $i (0 .. $#nucleotides) {
    my $mutations = ['A', 'G', 'C', 'T'];
    $new_triplets = $self->get_mutated_triplets($triplet, $mutations, $i, $new_triplets, $reverse);
  }
  return $new_triplets;
}

sub get_mutated_triplets {
  my $self = shift;
  my $triplet = shift;
  my $mutations = shift;
  my $position = shift;
  my $new_triplets = shift;
  my $reverse = shift;
  foreach my $mutation (@$mutations) {
    my $update_triplet = $triplet;
    if ($reverse) {
      my $reverse_mutation = $mutation;
      reverse_comp(\$reverse_mutation);
      substr($update_triplet, $position, 1, $reverse_mutation);
    } else {
      substr($update_triplet, $position, 1, $mutation);
    }
    $new_triplets->{$triplet}->{$position}->{$mutation} = $update_triplet;
  }
  return $new_triplets;
}

sub get_triplets {
  my $self = shift;
  my $translation_stable_id = shift;
  my $core_type = ($translation_stable_id =~ /^NP|XP/) ? 'otherfeatures' : 'core';
  my $cdba = $self->get_species_adaptor($core_type);
  my $translation_adaptor = $cdba->get_TranslationAdaptor or die "Failed to get translation adaptor";
  my $translation = $translation_adaptor->fetch_by_stable_id($translation_stable_id);
  my $slice_adaptor = $cdba->get_SliceAdaptor or die "Failed to get slice adaptor";

  my $transcript = $translation->transcript;
  my $chrom = $transcript->seq_region_name;
  my $start = $transcript->seq_region_start;
  my $end = $transcript->seq_region_end;
  my $strand = $transcript->seq_region_strand;
  my $slice = $slice_adaptor->fetch_by_region('toplevel', $chrom,  $start, $end);
  my $transcript_mapper = $transcript->get_TranscriptMapper();

  my $codonTable = Bio::Tools::CodonTable->new();
  my @all_triplets = ();
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
    my $entry = {
      coords => \@coords,
      aa_position => $i,
      chrom => $chrom,
      triplet_seq => $triplet,
    };
    my $aa = $codonTable->translate($triplet);
    if (!$aa) {
      $entry->{aa} = 'X';
    } else {
      $entry->{aa} = $aa;
      my $reverse = ($strand < 0);
      my $new_triplets = $self->mutate($triplet, $reverse);
      $entry->{new_triplets} = $new_triplets;
    }
    push @all_triplets, $entry;
  } 
  return \@all_triplets;
}


1;
