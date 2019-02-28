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

package Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::ImportMouse;

use warnings;
use strict;

use File::Path qw(make_path);
use File::stat;
use LWP::Simple;
use Data::Dumper; #TODO: remove if not needed
use Bio::EnsEMBL::Registry;
#use Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::Constants qw(GWAS NONE);
#use Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::BasePhenotypeAnnotation qw($variation_dba);
use Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::Constants qw(IMPC MGI NONE species);
use base ('Bio::EnsEMBL::Variation::Pipeline::PhenotypeAnnotation::MouseBasePhenotypeAnnotation');

my %source_info;
my $workdir;

my $debug;

sub fetch_input {
    my $self = shift;

    my $pipeline_dir = $self->required_param('pipeline_dir');
    my $run_type = $self->required_param('run_type');

    $debug        = $self->param('debug_mode');

    my $workdir = $pipeline_dir."/ImportMouse/";
    make_path($workdir);

    #get mouse coordinate file:
    my $coord_file = "MGI_MRK_Coord.rpt";
    my $impc_file_url = "http://www.informatics.jax.org/downloads/reports/MGI_MRK_Coord.rpt";
    print "Found file (".$workdir.$coord_file."), will skip new fetch\n" if -e $workdir."/".$coord_file;
    getstore($impc_file_url, $workdir."/".$coord_file) unless -e $workdir."/".$coord_file;

    unless ($run_type eq NONE) {
      my %import_species = &species;
      if($run_type eq IMPC){
        my @speciesList = map { {species => $_} } @{$import_species{'IMPC'}};
        foreach my $spec (@speciesList){
          $spec->{coord_file} = $workdir."/".$coord_file ;
          $spec->{pipeline_dir} = $workdir;
        }
        $self->param('output_ids', [ @speciesList ]);
        print "Setting up for IMPC import: ". join(", ",@{$import_species{'IMPC'}}). "\n" if $debug ;

      } elsif($run_type eq MGI){
        my @speciesList = map { {species => $_} } @{$import_species{'MGI'}};
        foreach my $spec (@speciesList){
          $spec->{coord_file} = $workdir."/".$coord_file;
          $spec->{pipeline_dir} = $workdir;
        }
        $self->param('output_ids', [ @speciesList ]);
        print "Setting up for MGI import: ". join(", ",@{$import_species{'MGI'}}). "\n" if $debug ;
      }
    }
}

sub write_output {
  my $self = shift;

  my $run_type = $self->param('run_type');
  unless ($run_type eq NONE) {
    if ($run_type eq IMPC){
      $self->dataflow_output_id($self->param('output_ids'), 2);
      print "Setting up for IMPC import: ".$self->param('output_ids')->{species}." species\n" if $self->param('debug_mode');
    } elsif ( $run_type eq MGI){
      $self->dataflow_output_id($self->param('output_ids'), 3);
      print "Setting up for MGI import: ".$self->param('output_ids')->{species}."\n" if $self->param('debug_mode');
    }
  }


}

1;
