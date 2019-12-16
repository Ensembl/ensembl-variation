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


=head1 ImportMouse

This module fetches the main coordinates file from http://www.informatics.jax.org.

=cut

package Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::ImportMouse;

use warnings;
use strict;

use File::Path qw(make_path);
use File::stat;
use LWP::Simple;

use Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::Constants qw(IMPC MGI MOUSE NONE SPECIES);
use base ('Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::MouseBasePhenotypeAnnotation');

sub fetch_input {
    my $self = shift;

    my $pipeline_dir = $self->required_param('pipeline_dir');
    my $run_type = $self->required_param('run_type');

    $self->debug($self->param('debug_mode'));

    my $workdir = $pipeline_dir."/ImportMouse/";
    make_path($workdir) or die "Failed to create $workdir $!\n";

    open(my $logFH, ">", $workdir."/".'log_import_out_ImportMouse_'.$run_type) || die ("Could not open file for writing $!\n");
    $self->logFH($logFH);

    #get mouse coordinate file:
    my $coord_file = "MGI_MRK_Coord.rpt";
    my $impc_file_url = "http://www.informatics.jax.org/downloads/reports/MGI_MRK_Coord.rpt";
    print $logFH "Found file (".$workdir.$coord_file."), will skip new fetch\n" if -e $workdir."/".$coord_file;
    getstore($impc_file_url, $workdir."/".$coord_file) unless -e $workdir."/".$coord_file;

    unless ($run_type eq NONE) {
      my %import_species = &SPECIES;
      if($run_type eq IMPC || $run_type eq MOUSE){
        my @speciesList = map { {species => $_} } @{$import_species{'IMPC'}};
        foreach my $spec (@speciesList){
          $spec->{coord_file} = $workdir."/".$coord_file ;
          $spec->{pipeline_dir} = $workdir;
          $spec->{source_name} = $run_type eq IMPC ? IMPC : 'IMPC_MGI' ;
        }
        $self->param('output_ids', [ @speciesList ]);
        print $logFH "Setting up for IMPC import: ". join(", ",@{$import_species{'IMPC'}}). "\n" if ($self->debug) ;

      } elsif($run_type eq MGI){
        my @speciesList = map { {species => $_} } @{$import_species{'MGI'}};
        foreach my $spec (@speciesList){
          $spec->{coord_file} = $workdir."/".$coord_file;
          $spec->{pipeline_dir} = $workdir;
        }
        $self->param('output_ids', [ @speciesList ]);
        print $logFH "Setting up for MGI import: ". join(", ",@{$import_species{'MGI'}}). "\n" if ($self->debug) ;
      }

      my @speciesNames = map { {species => $_->{species}} }  @{$self->param('output_ids')};
      $self->param('species_names', [@speciesNames]);
    }
}

sub write_output {
  my $self = shift;

  my $run_type = $self->param('run_type');
  unless ($run_type eq NONE) {
    if ($run_type eq IMPC || $run_type eq MOUSE){
      $self->print_logFH("Passing on import jobs (". scalar @{$self->param('output_ids')} .") for IMPC import \n") if ($self->debug);
      $self->dataflow_output_id($self->param('output_ids'), 2);
    } elsif ( $run_type eq MGI){
      $self->print_logFH("Passing on import jobs (". scalar @{$self->param('output_ids')} .") for MGI import \n") if ($self->debug);
      $self->dataflow_output_id($self->param('output_ids'), 3);
    }
    $self->print_logFH("Passing on check jobs (". scalar @{$self->param('output_ids')} .") for check_phenotypes \n") if ($self->debug);
    $self->dataflow_output_id($self->param('species_names'), 1);
  }
  close($self->logFH) if defined $self->logFH;

}

1;
