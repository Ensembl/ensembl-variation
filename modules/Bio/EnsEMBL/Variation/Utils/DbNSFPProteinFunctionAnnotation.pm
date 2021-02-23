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
with dbNSFP (revel, meta_lr and mutation_assessor) scores and predictions.

=cut


use strict;
use warnings;

package Bio::EnsEMBL::Variation::Utils::DbNSFPProteinFunctionAnnotation;
use Bio::EnsEMBL::Variation::Utils::BaseProteinFunctionAnnotation;
our @ISA = ('Bio::EnsEMBL::Variation::Utils::BaseProteinFunctionAnnotation');

my $REVEL_CUTOFF = 0.5;

=head2 new

  Example    :
  my $dbnsfp = Bio::EnsEMBL::Variation::Utils::DbNSFPProteinFunctionAnnotation->new(
    -species => 'Homo_sapiens',
    -annotation_file => 'dbNSFP3.5a_grch37.txt.gz',
    -assembly => 'GRCh37',
    -annotation_file_version => '3.5a',
    -pipeline_mode => 0,
    -debug_mode => 1,
  );

  Description: Constructor. Instantiates a new DbNSFPProteinFunctionAnnotation object.
  Returntype : DbNSFPProteinFunctionAnnotation
  Exceptions : throws on unsupported version
  Caller     : Bio::EnsEMBL::Variation::Pipeline::ProteinFunction::RunDbNSFP::run
  Status     : Stable
=cut
sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);

  if (! grep {$_ eq $self->annotation_file_version} ('3.5a', '4.0a', '4.1a')) {
    die "dbNSFP version " . $self->annotation_file_version . " is not supported.";
  }

  $self->analysis([qw/dbnsfp_revel dbnsfp_meta_lr dbnsfp_mutation_assessor/]);

  return $self;
}

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
#The rankscore cutoffs between "H" and "M", "M" and "L", and "L" and "N", are 0.9307, 0.52043 and 0.19675,

# dbNSFP has assembly and version specific headers
# If a new version is added, the new version first needs to tested and the column_names
# hash needs to be updated 
# the column names are changed for readability and consistency
my $column_names = {
  '3.5a' => {
    assembly_unspecific => {
      chr => '#chr',
      ref => 'ref',
      refcodon => 'refcodon',
      alt => 'alt',
      aaalt => 'aaalt',
      aaref => 'aaref',
      revel_score => 'REVEL_score',
      meta_lr_score => 'MetaLR_score',
      meta_lr_pred => 'MetaLR_pred',
      mutation_assessor_score => 'MutationAssessor_score_rankscore',
      mutation_assessor_pred => 'MutationAssessor_pred',
    },
    'assembly_specific' => {
      'GRCh37' => {
        pos => 'hg19_pos(1-based)'
      },
      'GRCh38' => {
        pos => 'pos(1-based)'
      },
    },
  },
  '4.0a' => {
    assembly_unspecific => {
      chr => '#chr',
      ref => 'ref',
      refcodon => 'refcodon',
      alt => 'alt',
      aaalt => 'aaalt',
      aaref => 'aaref',
      revel_score => 'REVEL_score',
      meta_lr_score => 'MetaLR_score',
      meta_lr_pred => 'MetaLR_pred',
      mutation_assessor_score => 'MutationAssessor_rankscore',
      mutation_assessor_pred => 'MutationAssessor_pred',
    },
    'assembly_specific' => {
      'GRCh37' => {
        pos => 'hg19_pos(1-based)'
      },
      'GRCh38' => {
        pos => 'pos(1-based)'
      },
    },
  },
  '4.1a' => {
    assembly_unspecific => {
      chr => '#chr',
      ref => 'ref',
      refcodon => 'refcodon',
      alt => 'alt',
      aaalt => 'aaalt',
      aaref => 'aaref',
      revel_score => 'REVEL_score',
      meta_lr_score => 'MetaLR_score',
      meta_lr_pred => 'MetaLR_pred',
      mutation_assessor_score => 'MutationAssessor_rankscore',
      mutation_assessor_pred => 'MutationAssessor_pred',
    },
    'assembly_specific' => {
      'GRCh37' => {
        pos => 'hg19_pos(1-based)'
      },
      'GRCh38' => {
        pos => 'pos(1-based)'
      },
    },
  }
};

sub load_predictions_for_triplets {
  my $self = shift;
  my $triplets = shift; 
  foreach my $entry (@$triplets) {
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
        my $data = $self->get_dbNSFP_row($line);
        my $chr = $data->{'chr'};
        my $pos = $data->{'pos'};
        my $ref = $data->{'ref'};
        my $refcodon = $data->{'refcodon'};
        my $alt = $data->{'alt'};
        my $aaalt = $data->{'aaalt'};
        my $aaref = $data->{'aaref'};
        next if ($alt eq $ref);
        my $nucleotide_position = ($self->reverse) ? $triplet_end - $pos : $pos - $triplet_start;
        my $mutated_triplet =  $new_triplets->{$triplet_seq}->{$nucleotide_position}->{$alt};
        my $mutated_aa = $self->codon_table->translate($mutated_triplet);
        next if ($aaalt ne $mutated_aa);
        $self->add_predictions($data, $i, $mutated_aa);
      }    
    }
  }
} 

sub add_predictions {
  my ($self, $data, $i, $mutated_aa) = @_;
  if ($data->{revel_score} ne '.') {
    my $prediction = ($data->{revel_score} >= $REVEL_CUTOFF) ? 'likely disease causing' : 'likely benign';
    $self->add_prediction($i, $mutated_aa, 'dbnsfp_revel', $data->{revel_score}, $prediction);
  }
  if ($data->{meta_lr_score} ne '.') {
    my $prediction = $predictions->{dbnsfp_meta_lr}->{$data->{meta_lr_pred}};
    $self->add_prediction($i, $mutated_aa, 'dbnsfp_meta_lr', $data->{meta_lr_score}, $prediction);
  }
  if ($data->{mutation_assessor_score} ne '.') {
    my $prediction;
    if ($self->annotation_file_version eq '3.5a') {
      $prediction = $predictions->{dbnsfp_mutation_assessor}->{$data->{mutation_assessor_pred}};  
    } elsif ($self->annotation_file_version eq '4.0a' || $self->annotation_file_version eq '4.1a') { 
      # In 4.0a the prediction is not always provided and we need to assign it based on the score thresholds     
      # The rankscore cutoffs between "H" and "M", "M" and "L", and "L" and "N", are 0.9307, 0.52043 and 0.19675,
      my $score = $data->{mutation_assessor_score}; 
      if ($score >= 0.9307) {
        $prediction = 'high';
      } elsif ($score >= 0.52043) {
        $prediction = 'medium'
      } elsif ($score >= 0.19675) {
        $prediction = 'low'
      } else {
        $prediction = 'neutral';
      }
    } else {
      die "dbNSFP version " . $self->annotation_file_version . " is not supported.";
    }
    $self->add_prediction($i, $mutated_aa, 'dbnsfp_mutation_assessor', $data->{mutation_assessor_score}, $prediction);
  }
}

=head2 get_dbNSFP_row

  Arg 1      : String $line from parser
  Description: - Join header column with row value
               - Use assembly and file version specific header
  Returntype : Hashref mapping header column to row value
  Exceptions : None
  Caller     : load_predictions_for_triplets()
  Status     :
=cut
sub get_dbNSFP_row {
  my $self = shift;
  my $line = shift;
  my @split = split /\t/, $line;
  my $header = $self->header;
  my $assembly = $self->assembly;
  my $dbnsfp_version = $self->annotation_file_version;
  my %raw_data = map {$header->[$_] => $split[$_]} (0..(scalar @{$header} - 1));
  my $data = {};
  if (!defined $column_names->{$dbnsfp_version}) {
    die "dbNSFP file column names are not specified for  version $dbnsfp_version";
  }
  my $assembly_unspecific = $column_names->{$dbnsfp_version}->{assembly_unspecific};
  foreach my $column_name (keys %{$assembly_unspecific}) {
    $data->{$column_name} = $raw_data{$assembly_unspecific->{$column_name}};
  }
  my $assembly_specific =  $column_names->{$dbnsfp_version}->{assembly_specific}->{$assembly}; 
  foreach my $column_name (keys %{$assembly_specific}) {
    $data->{$column_name} = $raw_data{$assembly_specific->{$column_name}};
  }
  return $data;
}

1;
