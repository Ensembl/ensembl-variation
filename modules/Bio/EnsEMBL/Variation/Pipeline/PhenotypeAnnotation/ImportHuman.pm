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


=head1 ImportHuman

This module triggers Human specific phenotype imports

=cut

package Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::ImportHuman;

use warnings;
use strict;

use Bio::EnsEMBL::Variation::Utils::SeqRegionUtils qw(update_seq_region_ids);
use Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::Constants qw(GWAS EGA ORPHANET MIMMORBID DDG2P CGC HUMAN HUMAN_VAR HUMAN_GENE NONE SPECIES);
use base ('Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::BasePhenotypeAnnotation');

# branch numbers from config file
my %source2branch = (
  GWAS      => 2,
  EGA       => 3,
  ORPHANET  => 4,
  MIMMORBID => 5,
  DDG2P     => 6,
  CGC       => 7,

  HUMAN     => 2,
  HUMAN_VAR => 2,
  HUMAN_GENE=> 4,
);

sub fetch_input {
    my $self = shift;

    my $pipeline_dir = $self->required_param('pipeline_dir');
    my $run_type = $self->required_param('run_type');

    open(my $pipelogFH, ">", $pipeline_dir."/".'log_import_debug_pipe_human') || die ("Failed to open file: $!\n");
    $self->pipelogFH($pipelogFH);

    unless ($run_type eq NONE) {
      my %import_species = SPECIES;
      #if HUMAN runtype, then select species by looking up MIMMORBID species
      # expectation is that MIMMORBID will always be only homo_sapiens
      die ("$run_type not defined in ImportHuman!\n") if (!$source2branch{$run_type});

      $self->param('output_ids',  [ map { {run_type => $run_type, species => $_} } @{$import_species{$run_type}} ]);
      $self->print_pipelogFH("Setting up for $run_type import: ". join(", ",@{$import_species{$run_type}}). "\n") if $self->param('debug_mode') ;
    }

    # if gene imports are performed then check first that the seq_region ids are in sync
    if ($run_type ne GWAS && $run_type ne EGA && $run_type ne HUMAN_VAR){
      # species parameter needs to be set for the core db, variation db adaptor fetch
      $self->param('species', 'homo_sapiens');
      update_seq_region_ids($self->core_db_adaptor, $self->variation_db_adaptor);
    }
}

sub write_output {
  my $self = shift;

  my $run_type = $self->param('run_type');
  unless ($run_type eq NONE) {
    if ($source2branch{$run_type}){
      $self->dataflow_output_id($self->param('output_ids'), $source2branch{$run_type});
      $self->print_pipelogFH( "Passing to $run_type import: ".scalar @{$self->param('output_ids')}." species\n") if $self->param('debug_mode');
    } else {
      $self->print_pipelogFH("Runtype $run_type not supported!\n");
      close($self->pipelogFH) if defined $self->pipelogFH;
      die("Unsupported $run_type\n");
    }

    $self->print_pipelogFH("Passing on check jobs (". scalar @{$self->param('output_ids')} .") for check_phenotypes \n") if $self->param('debug_mode');

    $self->dataflow_output_id($self->param('output_ids'), 8);
  }
  close($self->pipelogFH) if defined $self->pipelogFH;

}

1;
