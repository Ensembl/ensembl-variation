=head1 LICENSE

Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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
package Bio::EnsEMBL::Variation::Pipeline::DumpVEP::FinishDump;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Pipeline::DumpVEP::BaseVEP);

use File::Path;


sub param_defaults {
  return {
    'species_refseq' => 0,
    'merged' => 0,
  };
}

sub run {
  my $self = shift;
  my $refseq = $self->param('species_refseq') ? '_refseq' : '';
  
  $self->rm_dirs($refseq);
  $self->rm_dirs('_merged') if $refseq && $self->param('merged');
  
  return;
}

sub rm_dirs {
  my $self = shift;
  my $type = shift;
  $type ||= '';
  
  my $species  = $self->required_param('species');
  my $assembly = $self->required_param('assembly');
  my $version  = $self->required_param('ensembl_release');
  my $dir      = $self->required_param('pipeline_dir');

  rmtree(
    sprintf(
      '%s/%s%s/%s_%s',
      $dir,
      $species,
      $type,
      $version,
      $assembly
    )
  );
  
  # remove species dir if empty
  my $species_dir = "$dir/$species$type";
  
  opendir DIR, $species_dir or return;
  my @list = grep {!/^\.+$/} readdir DIR;
  closedir DIR;
  
  rmtree($species_dir) unless scalar @list;
  
  return;
}

1;
