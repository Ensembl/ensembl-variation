=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2025] EMBL-European Bioinformatics Institute

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

package Bio::EnsEMBL::Variation::Pipeline::InitRegulationEffect;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);

sub fetch_input {
  my $self = shift;
  my $species = $self->param_required('species');
  my $vdba = $self->get_adaptor($species, 'variation');
  my $vdbc = $vdba->dbc();
  $vdbc->do('TRUNCATE TABLE motif_feature_variation');
  $vdbc->do('ALTER TABLE motif_feature_variation DISABLE KEYS');
  $vdbc->do('TRUNCATE TABLE regulatory_feature_variation');
  $vdbc->do('ALTER TABLE regulatory_feature_variation DISABLE KEYS');
}

sub write_output {
  my $self = shift;
  return;
}
1;
