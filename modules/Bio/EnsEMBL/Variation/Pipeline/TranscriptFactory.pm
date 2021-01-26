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

package Bio::EnsEMBL::Variation::Pipeline::TranscriptFactory;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::TranscriptVariation;
use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(MAX_DISTANCE_FROM_TRANSCRIPT overlap);
use Bio::EnsEMBL::Variation::Utils::FastaSequence qw(setup_fasta);

use ImportUtils qw(load);
use FileHandle;
use Fcntl qw(:flock SEEK_END);
use Digest::MD5 qw(md5_hex);

use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);

my $DEBUG   = 0;

sub run {
  my $self = shift;

  my $gene_stable_id = $self->required_param('gene_stable_id'); 
  my $max_distance = $self->param('max_distance');

  # clear the registry here
  # this hopefully prevents any sequence caching issues
  # overhanging from previous jobs executed in the same hive process
  Bio::EnsEMBL::Registry->clear();

  my $core_dba = $self->get_species_adaptor('core');
  $core_dba->dbc->reconnect_when_lost(1);
  my $ga = $core_dba->get_GeneAdaptor;
  my $gene = $ga->fetch_by_stable_id($gene_stable_id)
      or die "failed to fetch gene for stable id: $gene_stable_id";


  my @transcript_ids_funnel = (); # e-hive speak: waits for fan jobs to finish
  my @transcript_ids_fan = (); # e-hive speak: fan transcript_ids (run for each transcript_id in parallel) 

  for my $transcript (@{ $gene->get_all_Transcripts }) { 
    push @transcript_ids_fan, {
      transcript_stable_id => $transcript->stable_id,
      max_distance => $max_distance,
      analysis => 'by_transcript',
    };
    push @transcript_ids_funnel, $transcript->stable_id,
  }

  $self->param('transcript_ids_funnel', {'transcript_ids_funnel' => \@transcript_ids_funnel, 'gene_stable_id' => $gene_stable_id});
  $self->param('transcript_ids_fan', \@transcript_ids_fan);

  return;
}

sub write_output {
  my $self = shift;
  $self->dataflow_output_id($self->param('transcript_ids_fan'), 2);
  $self->dataflow_output_id($self->param('transcript_ids_funnel'), 1);
}

1;
