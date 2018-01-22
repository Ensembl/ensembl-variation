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
package Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::BaseDataDumpsProcess;

use strict;

use base ('Bio::EnsEMBL::Hive::Process');

use Bio::EnsEMBL::Registry;

sub data_dir {
  my ($self,$species) = @_;
  my $data_dump_dir = $self->param('pipeline_dir');
  my $species_division = $self->param('species_division');
  # If division is defined append the pipeline_dir
  if ($species_division)
  {
    $data_dump_dir = $data_dump_dir."/".$species_division;
  }
  return $data_dump_dir;
}

sub get_all_species {
    my $self = shift;
    my $registry = 'Bio::EnsEMBL::Registry';
    $registry->load_all($self->param('registry_file'));
    my $vdbas = $registry->get_all_DBAdaptors(-group => 'variation');
    my $species = {};
    foreach my $vdba (@$vdbas) {
        my $species_name = $vdba->species();
        $species->{$species_name} = 1;
    }
    return $species;
}

sub get_species_adaptor {
    my ($self, $species, $group) = @_;
    return $self->get_adaptor($species, $group);
}

sub get_adaptor {
    my ($self, $species, $group) = @_;
    my $dba;
    eval {
        $dba = Bio::EnsEMBL::Registry->get_DBAdaptor($species, $group);
    };
    unless (defined $dba) {
        $self->_load_registry();
        $dba = Bio::EnsEMBL::Registry->get_DBAdaptor($species, $group);
    }
    unless (defined $dba) {
        die "Failed to a get DBA for $species and group $group";
    }
    return $dba;
}

sub _load_registry {
    my ($self) = @_;
    my $reg_file = $self->param('registry_file');
    Bio::EnsEMBL::Registry->load_all($reg_file, 0, 1);
    return;
}


1;

