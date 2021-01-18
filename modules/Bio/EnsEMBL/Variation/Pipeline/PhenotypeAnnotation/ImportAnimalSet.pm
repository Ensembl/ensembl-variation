=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2021] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     https://www.apache.org/licenses/LICENSE-2.0

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


=head1 ImportAnimalSet

This module triggers a set of specific animal imports

=cut

package Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::ImportAnimalSet;

use warnings;
use strict;

use Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::Constants qw(ANIMALQTL OMIA ANIMALSET NONE SPECIES);
use base ('Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::BasePhenotypeAnnotation');

my %source2branch = (
  OMIA      => 2,
  ANIMALQTL => 3,
);

my $type;

sub fetch_input {
    my $self = shift;

    my $pipeline_dir = $self->required_param('pipeline_dir');
    my $run_type = $self->required_param('run_type');

    open(my $pipelogFH, ">", $pipeline_dir."/".'log_import_debug_pipe_animalSet') || die ("Failed to open file: $!\n");
    $self->pipelogFH($pipelogFH);

    unless ($run_type eq NONE) {
      my %import_species = SPECIES;
      $type = ($run_type eq ANIMALSET) ? 'OMIA' : $run_type;
      $self->param('output_ids',  [ map { {run_type => $run_type, species => $_} } @{$import_species{$type}} ]);
      $self->print_pipelogFH("Setting up for $run_type import: ". join(", ",@{$import_species{$type}}). "\n") if $self->param('debug_mode') ;
    }
}

sub write_output {
  my $self = shift;

  my $run_type = $self->param('run_type');
  unless ($run_type eq NONE) {
    if ($source2branch{$type}){
      $self->dataflow_output_id($self->param('output_ids'), $source2branch{$type});
      $self->print_pipelogFH( "Passing to $type import: ".scalar @{$self->param('output_ids')}." species\n") if $self->param('debug_mode');
    } else {
      $self->print_pipelogFH("Runtype $type not supproted!\n");
    }
  }
  close($self->logFH) if defined $self->logFH;

}

1;
