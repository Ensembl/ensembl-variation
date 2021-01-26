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

package Bio::EnsEMBL::Variation::Pipeline::DumpVariationGeneName;

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

  my $gene_id = $self->required_param('gene_stable_id'); 

  # clear the registry here
  # this hopefully prevents any sequence caching issues
  # overhanging from previous jobs executed in the same hive process
  Bio::EnsEMBL::Registry->clear();

  my $core_dba = $self->get_species_adaptor('core');
  $core_dba->dbc->reconnect_when_lost(1);

  my $var_dba = $self->get_species_adaptor('variation');
  $var_dba->dbc->reconnect_when_lost(1);
  
  my $ga = $core_dba->get_GeneAdaptor;
  my $sa = $core_dba->get_SliceAdaptor;


  my $gene = $ga->fetch_by_stable_id($gene_id) 
    or die "failed to fetch gene for stable id: $gene_id";

  my $gene_name = $gene->display_xref->display_id if defined $gene->display_xref();
  if (!$gene_name) {
    $self->warning("No display xref for $gene_id");
  }

  my $max_distance = $self->get_max_distance;

  my $slice = $sa->fetch_by_gene_stable_id(
    $gene->stable_id, 
    $max_distance
  ) or die "failed to get slice around gene: ".$gene->stable_id;
  
  # call seq here to help cache
  $slice->seq();

  $gene = $gene->transfer($slice);

  my $vfa = $var_dba->get_VariationFeatureAdaptor();

  my @vfs = (
    @{ $vfa->fetch_all_by_Slice_SO_terms($slice) },
    @{ $vfa->fetch_all_somatic_by_Slice_SO_terms($slice) }
  );

  my $web_index_files_dir = $self->get_files_dir($gene_id, 'web_index');

  $self->param('gene_stable_id', {'gene_stable_id' => $gene->stable_id, 'max_distance' => $max_distance});

  if (($gene->length > 1e6) ||
      (scalar(@vfs) > 500_000) ||
      (scalar(@{$gene->get_all_Transcripts()}) >= 50)) {
    $self->param('is_big_gene', 1);
  } 

  return;
}

sub write_output {
  my $self = shift;
  my $input = $self->param('gene_stable_id');
  if ($self->param('is_big_gene')) {
    $self->dataflow_output_id($self->param('gene_stable_id'), 2);
  } else {
    $self->dataflow_output_id($self->param('gene_stable_id'), 3);
  }
}

sub get_max_distance {
  my $self = shift;
  my $opt_max_distance = $self->param('max_distance');
  my $max_distance;

  if(defined($opt_max_distance)) {
    $max_distance = $opt_max_distance;
    $Bio::EnsEMBL::Variation::Utils::VariationEffect::UPSTREAM_DISTANCE = $opt_max_distance;
    $Bio::EnsEMBL::Variation::Utils::VariationEffect::DOWNSTREAM_DISTANCE = $opt_max_distance;
  }
  else {
    $max_distance = MAX_DISTANCE_FROM_TRANSCRIPT;
  }
  return $max_distance;
}

1;
