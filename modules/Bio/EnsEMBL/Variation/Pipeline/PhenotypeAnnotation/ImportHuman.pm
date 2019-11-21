=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2019] EMBL-European Bioinformatics Institute

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

This module triggers Human specific imports

=cut

package Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::ImportHuman;

use warnings;
use strict;

use File::Path qw(make_path);
use File::stat;
use LWP::Simple;

use Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::Constants qw(GWAS EGA Orphanet MIMmorbid DDG2P CGC Human NONE species);
use base ('Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::BasePhenotypeAnnotation');

sub fetch_input {
    my $self = shift;

    my $pipeline_dir = $self->required_param('pipeline_dir');
    my $run_type = $self->required_param('run_type');

    open (my $pipelogFH, ">", $pipeline_dir."/".'log_import_debug_pipe_human') || die ("Failed to open file: $!\n");
    $self->pipelogFH($pipelogFH);

    unless ($run_type eq NONE) {
      my %import_species = &species;
      $run_type = GWAS if $run_type eq Human; #GWAS analysis will trigger all the subsequent ones
      $self->param('output_ids',  [ map { {species => $_} } @{$import_species{$run_type}} ]);
      $self->print_pipelogFH("Setting up for $run_type import: ". join(", ",@{$import_species{$run_type}}). "\n") if $self->param('debug_mode') ;
    }
}

sub write_output {
  my $self = shift;

  my $run_type = $self->param('run_type');
  unless ($run_type eq NONE) {
    if ( $run_type eq GWAS || $run_type eq Human){
      $self->dataflow_output_id($self->param('output_ids'), 2);
      $self->print_pipelogFH( "Passing to NHGRI-EBI GWAS import: ".scalar @{$self->param('output_ids')}." species\n") if $self->param('debug_mode');
    } elsif ( $run_type eq EGA){
      $self->dataflow_output_id($self->param('output_ids'), 3);
      $self->print_pipelogFH( "Passing to EGA import: ".scalar @{$self->param('output_ids')}." species\n") if $self->param('debug_mode');
    } elsif ( $run_type eq Orphanet){
      $self->dataflow_output_id($self->param('output_ids'), 4);
      $self->print_pipelogFH( "Passing to Orphanet import: ".scalar @{$self->param('output_ids')}." species\n") if $self->param('debug_mode');
    } elsif ( $run_type eq MIMmorbid){
      $self->dataflow_output_id($self->param('output_ids'), 5);
      $self->print_pipelogFH( "Passing to MIMmorbid import: ".scalar @{$self->param('output_ids')}." species\n") if $self->param('debug_mode');
    } elsif ( $run_type eq DDG2P){
      $self->dataflow_output_id($self->param('output_ids'), 6);
      $self->print_pipelogFH( "Passing to DDG2P import: ".scalar @{$self->param('output_ids')}." species\n") if $self->param('debug_mode');
    } elsif ( $run_type eq CGC){
      $self->dataflow_output_id($self->param('output_ids'), 7);
      $self->print_pipelogFH( "Passing to CGC import: ".scalar @{$self->param('output_ids')}." species\n") if $self->param('debug_mode');
    }

    $self->print_pipelogFH("Passing on check jobs (". scalar @{$self->param('output_ids')} .") for check_phenotypes \n") if $self->param('debug_mode');
    $self->dataflow_output_id($self->param('output_ids'), 1);
  }
  close($self->logFH) if defined $self->logFH;

}

1;
