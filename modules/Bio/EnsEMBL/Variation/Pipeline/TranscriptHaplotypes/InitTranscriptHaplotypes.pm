=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

package Bio::EnsEMBL::Variation::Pipeline::TranscriptHaplotypes::InitTranscriptHaplotypes;

use File::Path qw(rmtree);

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);

my $DEBUG = 0;

sub fetch_input {
  my $self = shift;

  my $core_dba = $self->get_species_adaptor('core');

  my $ta = $core_dba->get_TranscriptAdaptor or die "Failed to get transcript adaptor";

  my (@transcript_output_ids, @big_transcript_output_ids);
  my $tr_count = 0;

  # my @transcripts = ($ta->fetch_by_stable_id('ENST00000304748'));
  my @transcripts = @{ $ta->fetch_all_by_biotype('protein_coding') };

  for my $tr (@transcripts) {
    next unless $tr->translation;

    $tr_count++;

    if($tr->translation->length > 2500) {
      push @big_transcript_output_ids, {
        transcript_stable_id => $tr->stable_id,
      }
    }

    else {
      push @transcript_output_ids, {
        transcript_stable_id  => $tr->stable_id,
      };
    }
    
    if ($DEBUG) {
      last if $tr_count >= 100;
    }
  }

  # push @transcript_output_ids, { transcript_stable_id => 'ENST00000595042' };

  if (@transcript_output_ids) {
    $self->param('transcript_output_ids', \@transcript_output_ids);
  }

  if (@big_transcript_output_ids) {
    $self->param('big_transcript_output_ids', \@big_transcript_output_ids);
  }

  my $dir = $self->required_param('pipeline_dir')."/haplotypes";
  rmtree($dir) if -d $dir;
  mkdir($dir) or die "ERROR: Could not make directory $dir\n";
}

sub write_output {
  my $self = shift;

  if (my $transcript_output_ids = $self->param('transcript_output_ids')) {
    $self->dataflow_output_id($transcript_output_ids, 1);
  }

  if (my $big_transcript_output_ids = $self->param('big_transcript_output_ids')) {
    $self->dataflow_output_id($big_transcript_output_ids, 2);
  }

  $self->dataflow_output_id([{}], 3);

  return;
}

1;
