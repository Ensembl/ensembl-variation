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

package Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::InitPhenotypeAnnotation;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::Constants qw(RGD AnimalQTL ZFIN GWAS OMIA EGA Orphanet MIMmorbid DDG2P CGC IMPC MGI NONE species);
use base qw(Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::BasePhenotypeAnnotation);

my $logPipeFH;
sub fetch_input {
  my $self = shift;

  my $pipeline_dir = $self->required_param('pipeline_dir');
  my $run_type = $self->required_param('run_import_type');

  my $debug = $self->param('debug_mode');

  open ($logPipeFH, ">", $pipeline_dir."/".'log_import_debug_pipe');
  unless ($run_type eq NONE) {
    my %import_species = &species;
    if ($run_type eq RGD){
      $self->param('output_ids',  [ map { {species => $_} } @{$import_species{'RGD'}} ]);
      print $logPipeFH "Setting up for RGD import: ". join(", ",@{$import_species{'RGD'}}). "\n" if $debug ;

    } elsif($run_type eq AnimalQTL){
      $self->param('output_ids',  [ map { {species => $_} } @{$import_species{'AnimalQTL'}} ]);
      print $logPipeFH "Setting up for AnimalQTL import: ". join(", ",@{$import_species{'AnimalQTL'}}). "\n" if $debug ;

    } elsif($run_type eq ZFIN){
      $self->param('output_ids',  [ map { {species => $_} } @{$import_species{'ZFIN'}} ]);
      print $logPipeFH "Setting up for ZFIN import: ". join(", ",@{$import_species{'ZFIN'}}). "\n" if $debug ;

    } elsif($run_type eq GWAS){
      $self->param('output_ids',  [ map { {species => $_} } @{$import_species{'GWAS'}} ]);
      print $logPipeFH "Setting up for NHGRI-EBI GWAS Catalog import: ". join(", ",@{$import_species{'GWAS'}}). "\n" if $debug ;

    } elsif($run_type eq OMIA){
      $self->param('output_ids',  [ map { {species => $_} } @{$import_species{'OMIA'}} ]);
      print $logPipeFH "Setting up for OMIA import: ". join(", ",@{$import_species{'OMIA'}}). "\n" if $debug ;

    } elsif($run_type eq EGA){
      $self->param('output_ids',  [ map { {species => $_} } @{$import_species{'EGA'}} ]);
      print $logPipeFH "Setting up for EGA import: ". join(", ",@{$import_species{'EGA'}}). "\n" if $debug ;

    } elsif($run_type eq Orphanet){
      $self->param('output_ids',  [ map { {species => $_} } @{$import_species{'Orphanet'}} ]);
      print $logPipeFH "Setting up for Orphanet import: ". join(", ",@{$import_species{'Orphanet'}}). "\n" if $debug ;

    } elsif($run_type eq MIMmorbid){
      $self->param('output_ids',  [ map { {species => $_} } @{$import_species{'MIMmorbid'}} ]);
      print $logPipeFH "Setting up for MIMmorbid import: ". join(", ",@{$import_species{'MIMmorbid'}}). "\n" if $debug ;

    } elsif($run_type eq DDG2P){
      $self->param('output_ids',  [ map { {species => $_} } @{$import_species{'DDG2P'}} ]);
      print $logPipeFH "Setting up for DDG2P import: ". join(", ",@{$import_species{'DDG2P'}}). "\n" if $debug ;

    } elsif($run_type eq CGC){
      $self->param('output_ids',  [ map { {species => $_} } @{$import_species{'CGC'}} ]);
      print $logPipeFH "Setting up for CancerGeneConsensus import: ". join(", ",@{$import_species{'CGC'}}). "\n" if $debug ;

    } elsif($run_type eq IMPC || $run_type eq MGI ){
      $self->param('output_ids', [{run_type => $run_type}]);
    } else {
      print $logPipeFH "WARNING: No valid run_import_type specified: $run_type\n" if $debug ;
    }
  }
}

sub write_output {
  my $self = shift;

  my $run_type = $self->param('run_import_type');

  unless ($run_type eq NONE) {
    if ($run_type eq RGD){
      $self->dataflow_output_id($self->param('output_ids'), 2);
      print $logPipeFH "Passing to RGD import: ".scalar @{$self->param('output_ids')}." species\n" if $self->param('debug_mode');
    } elsif ( $run_type eq AnimalQTL){
      $self->dataflow_output_id($self->param('output_ids'), 3);
      print $logPipeFH "Passing to AnimalQTL import: ".scalar @{$self->param('output_ids')}." species\n" if $self->param('debug_mode');
    } elsif ( $run_type eq ZFIN){
      $self->dataflow_output_id($self->param('output_ids'), 4);
      print $logPipeFH "Passing to ZFIN import: ".scalar @{$self->param('output_ids')}." species\n" if $self->param('debug_mode');
    } elsif ( $run_type eq GWAS){
      $self->dataflow_output_id($self->param('output_ids'), 5);
      print $logPipeFH "Passing to NHGRI-EBI GWAS import: ".scalar @{$self->param('output_ids')}." species\n" if $self->param('debug_mode');
    } elsif ( $run_type eq OMIA){
      $self->dataflow_output_id($self->param('output_ids'), 6);
      print $logPipeFH "Passing to OMIA import: ".scalar @{$self->param('output_ids')}." species\n" if $self->param('debug_mode');
    } elsif ( $run_type eq EGA){
      $self->dataflow_output_id($self->param('output_ids'), 7);
      print $logPipeFH "Passing to EGA import: ".scalar @{$self->param('output_ids')}." species\n" if $self->param('debug_mode');
    } elsif ( $run_type eq Orphanet){
      $self->dataflow_output_id($self->param('output_ids'), 8);
      print $logPipeFH "Passing to Orphanet import: ".scalar @{$self->param('output_ids')}." species\n" if $self->param('debug_mode');
    } elsif ( $run_type eq MIMmorbid){
      $self->dataflow_output_id($self->param('output_ids'), 9);
      print $logPipeFH "Passing to MIMmorbid import: ".scalar @{$self->param('output_ids')}." species\n" if $self->param('debug_mode');
    } elsif ( $run_type eq DDG2P){
      $self->dataflow_output_id($self->param('output_ids'), 10);
      print $logPipeFH "Passing to DDG2P import: ".scalar @{$self->param('output_ids')}." species\n" if $self->param('debug_mode');
    } elsif ( $run_type eq CGC){
      $self->dataflow_output_id($self->param('output_ids'), 11);
      print $logPipeFH "Passing to CancerGeneConsensus import: ".scalar @{$self->param('output_ids')}." species\n" if $self->param('debug_mode');
    } elsif ( $run_type eq IMPC || $run_type eq MGI){
      $self->dataflow_output_id($self->param('output_ids'), 12);
      print $logPipeFH "Passing to $run_type import \n" if $self->param('debug_mode');
    }
  }
  close($logPipeFH);
}

1;
