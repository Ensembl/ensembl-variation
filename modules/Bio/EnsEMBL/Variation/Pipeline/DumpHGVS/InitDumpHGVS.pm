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

=head1 NAME

Bio::EnsEMBL::Variation::Pipeline::DumpHGVS::InitDumpHGVS

=head1 DESCRIPTION

Initialise the HGVS dump

=cut

package Bio::EnsEMBL::Variation::Pipeline::DumpHGVS::InitDumpHGVS;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);

use FileHandle;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::Slice qw(split_Slices);
use File::Path qw(make_path);


sub fetch_input {
  my $self = shift;

  my $region_size = $self->param_required('region_size');
  my $region_overlap = $self->param_required('region_overlap');

  my @regions;

  my $core_dba = $self->get_species_adaptor('core');
  my $sa = $core_dba->get_SliceAdaptor();

  # Get all seq regions
  my $slices = $sa->fetch_all('chromosome');

  $self->warning('Number of chromosomes print (' . scalar(@$slices) . ') to do');

  my $sub_slices = split_Slices($slices, $region_size, $region_overlap);
  $self->warning('Number of subslices (' . scalar(@$sub_slices) . ') to do');

  for my $slice (@$sub_slices) {
    push @regions , { region => join(':', $slice->seq_region_name(),
                                join('-', $slice->start(), $slice->end()))};
  }
  $self->param('regions', \@regions);
}

sub run {
  my $self = shift;
}


sub write_output {
  my $self = shift;
  my $regions = $self->param('regions');
  $self->warning(scalar @{$regions} . ' jobs to do');
  $self->dataflow_output_id($regions, 2);
}

1;
