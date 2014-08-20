=head1 LICENSE

Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

=head1 NAME

Bio::EnsEMBL::Variation::ProteinHaplotype

=head1 SYNOPSIS

  use Bio::EnsEMBL::Variation::ProteinHaplotype;

=head1 DESCRIPTION

A class for representing a transcript's CDS sequence modified by individual
genotypes.

=cut

package Bio::EnsEMBL::Variation::ProteinHaplotype;

use strict;
use warnings;
use Bio::EnsEMBL::Variation::TranscriptHaplotype;

use base qw(Bio::EnsEMBL::Variation::TranscriptHaplotype);

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  
  my %args = @_;
  
  my $self = $class->SUPER::new(%args);
  bless($self, $class);
  
  return $self;
}

sub get_all_CDSHaplotypes {
  return $_[0]->get_other_Haplotypes;
}

sub get_all_diffs {
  my $self = shift;
  
  if(!defined($self->{diffs})) {
    my $missense_preds = $self->container->_get_missense_predictions;
    
    my @diffs;
    
    foreach my $raw_diff(@{$self->_get_raw_diffs}) {
      my $diff = {
        diff => $raw_diff,
      };
      
      for(qw(sift_score sift_prediction polyphen_score polyphen_prediction)) {
        $diff->{$_} = $missense_preds->{$raw_diff}->{$_} if defined($missense_preds->{$raw_diff}) && defined($missense_preds->{$raw_diff}->{$_});
      }
      
      push @diffs, $diff;
    }
    
    $self->{diffs} = \@diffs;
  }
  
  return $self->{diffs};
}

1;