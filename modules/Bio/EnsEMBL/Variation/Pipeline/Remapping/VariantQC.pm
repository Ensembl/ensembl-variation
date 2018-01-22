#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.




=head1 CONTACT

Please email comments or questions to the public Ensembl
developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

Questions may also be sent to the Ensembl help desk at
<helpdesk.org>.

=cut

package Bio::EnsEMBL::Variation::Pipeline::Remapping::VariantQC;

use base ('Bio::EnsEMBL::Hive::Process');

use FileHandle;

=begin
sub run {
  my $self = shift;
  my $file_number = $self->param('file_number');
  # failure reasons:
  # None of the variant alleles match the reference allele 
  # Variation can not be re-mapped to the current assembly
  # Variation maps to more than one genomic location
  my $fh = FileHandle->new("", 'r');
  my $failed_variation_ids = {};

  while (<$fh>)  {
    chomp;
    my $data = $self->read_line($_);
    my $variation_id = $data->{variation_id};
    my $seq_region_id = $data->{seq_region_id};
    my $seq_region_start = $data->{seq_region_start};
    my $seq_region_end = $data->{seq_region_end};
    my $map_weight = $data->{map_weight};
    my $allele_string = $data->{allele_string};
    if ($map_weight > 1) {
      $description = 'Variation maps to more than one genomic location';
    }
    my $allele_maps_reference = 0
    if (!$allele_maps_reference) {
      $description = 'None of the variant alleles match the reference allele';
    }
  } 
}

sub set_failed_description {
  my $self = shift;
  # failure reasons:
  # None of the variant alleles match the reference allele 
  # Variation can not be re-mapped to the current assembly
  # Variation maps to more than one genomic location
}
=end
=cut


1;
