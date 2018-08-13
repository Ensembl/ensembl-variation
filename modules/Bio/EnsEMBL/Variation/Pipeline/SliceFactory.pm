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

package Bio::EnsEMBL::Variation::Pipeline::SliceFactory;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);

sub fetch_input {
    my $self = shift;
    my $species = $self->param_required('species');

    my $cdba = $self->get_adaptor($species, 'core');

    my $sa = $cdba->get_SliceAdaptor or die 'Failed to get SliceAdaptor';

    my $slices = $sa->fetch_all('toplevel', undef, 0, 1);
    my @slice_names = ();

   if ($self->param('debug')) {
      my $slice = $sa->fetch_by_region('chromosome', 12);
      push @slice_names, {slice_name => $slice->seq_region_name };
    } else {
      foreach my $slice (@$slices) {
        push @slice_names, { slice_name => $slice->seq_region_name}; 
      } 
    }

    $self->dataflow_output_id(\@slice_names, 2);
}


sub write_output {
    my $self = shift;
    return;
}
1;
