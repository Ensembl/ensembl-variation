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

package Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::InitPhenotypeAnnotation;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::Constants qw(RGD ANIMALQTL ZFIN GWAS OMIA EGA ORPHANET MIMMORBID DDG2P CGC IMPC MGI MOUSE HUMAN ANIMALSET NONE SPECIES GROUP_RUN_TYPES SOURCES_IN_RUN_TYPES);

use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);

my $logPipeFH;

my %source2branch = (
  HUMAN     => 2,
  HUMAN_VAR => 2,
  HUMAN_GENE=> 2,
  GWAS      => 2,
  EGA       => 2,
  ORPHANET  => 2,
  MIMMORBID => 2,
  DDG2P     => 2,
  CGC       => 2,

  MOUSE     => 3,
  IMPC      => 3,
  MGI       => 3,

  ANIMALSET => 4,
  OMIA      => 4,
  ANIMALQTL => 4,

  RGD       => 5,
  ZFIN      => 6
);


sub fetch_input {
  my $self = shift;

  my $pipeline_dir = $self->required_param('pipeline_dir');
  my $run_type = $self->required_param('run_type');
  $run_type =uc($run_type);

  my $debug = $self->param('debug_mode');

  open($logPipeFH, ">", $pipeline_dir."/".'log_import_debug_pipe') || die ("Failed to open file: $!\n");

  # fetch species to be annotated
  unless ($run_type eq NONE) {
    my %import_species = SPECIES;
    my %group_runs = GROUP_RUN_TYPES;
    my %source_runs = SOURCES_IN_RUN_TYPES;

    #MOUSE, HUMAN and specific AnimalSet have separate init analysis
    if ($group_runs{$run_type} || $source_runs{$run_type}){
      $self->param('output_ids', [{run_type => $run_type}]);
    } elsif ($import_species{$run_type}) {
      $self->param('output_ids',  [ map { {run_type=> $run_type, species => $_} } @{$import_species{$run_type}} ]);
      print $logPipeFH "Setting up for $run_type import: ". join(", ",@{$import_species{$run_type}}). "\n" if $debug ;
    } else {
      print $logPipeFH "WARNING: No valid run_import_type specified: $run_type\n" if $debug ;
    }
  }
}

sub write_output {
  my $self = shift;

  my $run_type = $self->param('run_type');
  $run_type =uc($run_type);

  unless ($run_type eq NONE) {
    if ($source2branch{$run_type}){
      $self->dataflow_output_id($self->param('output_ids'), $source2branch{$run_type});
      print $logPipeFH "Passing to $run_type import: ".scalar @{$self->param('output_ids')}." species\n" if $self->param('debug_mode');
    } else {
      print $logPipeFH "Run type $run_type not supported!\n";
    }
  }
  close($logPipeFH);
}

1;
