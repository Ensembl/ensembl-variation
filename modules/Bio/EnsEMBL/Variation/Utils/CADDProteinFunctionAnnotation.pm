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

=head1 DESCRIPTION

Module is used in protein function prediction pipeline for
annotating all possible amino acid substitutions in a translation
with CADD scores and predictions. 

=cut


use strict;
use warnings;

package Bio::EnsEMBL::Variation::Utils::CADDProteinFunctionAnnotation;

use Bio::EnsEMBL::Variation::Utils::BaseProteinFunctionAnnotation;
our @ISA = ('Bio::EnsEMBL::Variation::Utils::BaseProteinFunctionAnnotation');

my $CADD_CUTOFF = 30;

=head2 new
  Example    :
  my $cadd = Bio::EnsEMBL::Variation::Utils::CADDProteinFunctionAnnotation->new(
    -species => 'Homo_sapiens',
    -annotation_file => 'cadd_v1.3_grch37.txt.gz',
    -assembly => 'GRCh37',
    -annotation_file_version => 'v1.3',
    -pipeline_mode => 0,
    -debug_mode => 1,
  );

  Description: Constructor. Instantiates a new CADDProteinFunctionAnnotation object.
  Returntype : CADDProteinFunctionAnnotation
  Exceptions : throws on unsupported version
  Caller     : Bio::EnsEMBL::Variation::Pipeline::ProteinFunction::RunCADD::run
  Status     : Stable
=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);

  if (! grep {$_ eq $self->annotation_file_version} ('v1.3', 'v1.4', 'v1.5', 'v1.6')) {
    die "CADD version " . $self->annotation_file_version . " is not supported.";
  }

  $self->analysis([qw/cadd/]);

  return $self;
}

sub load_predictions_for_triplets {
  my $self = shift;
  my $all_triplets = shift; 
  foreach my $entry (@$all_triplets) {
    my $aa = $entry->{aa};
    $self->amino_acids($aa);
    next if $aa eq 'X';
    my @coords = @{$entry->{coords}};
    my $chrom = $entry->{chrom};
    my $triplet_seq = $entry->{triplet_seq};
    my $i = $entry->{aa_position};
    my $new_triplets = $entry->{new_triplets};
    foreach my $coord (@coords) {
      my $triplet_start = $coord->[0];
      my $triplet_end = $coord->[1];
      my $iter = $self->get_tabix_iterator($chrom, $triplet_start, $triplet_end);
      next if (!defined $iter);
      while (my $line = $iter->next) {
        my $data = $self->get_CADD_row($line);
        my $chr = $data->{'#Chr'};
        my $pos = $data->{'Pos'};
        my $ref = $data->{'Ref'};
        my $alt = $data->{'Alt'};
        my $cadd_phred = $data->{'PHRED'};
        next if ($alt eq $ref);
        my $nucleotide_position = ($self->reverse) ? $triplet_end - $pos : $pos - $triplet_start;
        my $mutated_triplet =  $new_triplets->{$triplet_seq}->{$nucleotide_position}->{$alt};
        my $mutated_aa = $self->codon_table->translate($mutated_triplet);
        $self->add_predictions($data, $i, $mutated_aa);
      }
    }
  }
} 

sub add_predictions {
  my ($self, $data, $i, $mutated_aa) = @_;
  if ($data->{'PHRED'} ne '.') {
    my $prediction = ($data->{'PHRED'} >= $CADD_CUTOFF) ? 'likely deleterious' : 'likely benign';
    $self->add_prediction($i, $mutated_aa, 'cadd', $data->{'PHRED'}, $prediction);
  }
}

=head2 get_CADD_row
  Arg 1      : String $line from parser
  Description: Join header column with row value
  Returntype : Hashref mapping header column to row value
  Exceptions : None
  Caller     : load_predictions_for_triplets() 
  Status     : 
=cut
sub get_CADD_row {
  my $self = shift;
  my $line = shift;
  my @split = split /\t/, $line;
  my $header = $self->header;
  my %data = map {$header->[$_] => $split[$_]} (0..(scalar @{$header} - 1));
  return \%data;
}

1;
