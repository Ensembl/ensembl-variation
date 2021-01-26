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

Bio::EnsEMBL::Variation::Pipeline::GetCAR::InitGetCAR

=head1 DESCRIPTION

Initialise GetCAR pipeline

=cut

package Bio::EnsEMBL::Variation::Pipeline::GetCAR::InitGetCAR;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);

sub fetch_input {
  my $self = shift;

  my $hgvs_dir = $self->param_required('hgvs_dir');
  my $car_lu_dir = $self->param_required('car_lu_dir');

  if (! -e $hgvs_dir) {
    die("hgvs_dir ($hgvs_dir) does not exist");
  }
  if (! -e $car_lu_dir) {
     die("car_lu_dir ($car_lu_dir) does not exist");
  }
  my @seq_names = (1..22, 'X','Y', 'MT');
  @seq_names = map {{seq_name => $_}} @seq_names;

  $self->param('seq_names', \@seq_names);
}

sub run {
  my $self = shift;
}


sub write_output {
  my $self = shift;
  my $seq_names = $self->param('seq_names');
  $self->warning(scalar @{$seq_names} . ' jobs to do');
  $self->dataflow_output_id($seq_names, 2);
}

1;
