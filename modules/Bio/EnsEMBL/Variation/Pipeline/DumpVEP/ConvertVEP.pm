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
package Bio::EnsEMBL::Variation::Pipeline::DumpVEP::ConvertVEP;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Pipeline::DumpVEP::BaseVEP);


sub param_defaults {
  return {
    'merged' => 0,
    'species_refseq' => 0,
  };
}

sub run {
  my $self = shift;
  
  my $type = $self->param('species_refseq') ? 'refseq' : '';
  
  $self->convert($type);
  
  # convert the merged cache too if this is the RefSeq one
  $self->convert('merged', '_tabixconverted') if $type;
}

sub convert {
  my $self = shift;
  my $type = shift;
  my $mod  = shift;
  $type ||= '';
  
  # basic params
  my $debug    = $self->param('debug');
  my $species  = $self->required_param('species');
  my $assembly = $self->required_param('assembly');
  my $version  = $self->required_param('ensembl_release');
  my $dir      = $self->required_param('pipeline_dir');
  
  $species .= '_'.$type if $type;
  
  # cmnd line
  my $perl    = $self->required_param('perl_command');
  my $con_cmd = $self->required_param('ensembl_cvs_root_dir').'/ensembl-tools/scripts/variant_effect_predictor/convert_cache.pl';
  
  # construct command
  my $cmd = sprintf(
    '%s %s --species %s --version %s_%s --dir %s --remove',
    $perl,
    $con_cmd,
    
    $species,
    $version,
    $assembly,
    $dir,
  );
  
  my $finished = 0;
  
  if($debug) {
    print STDERR "$cmd\n";
  }
  else {
    open CMD, "$cmd 2>&1 |" or die "ERROR: Failed to run command $cmd";
    my @buffer;
    while(<CMD>) {
      $finished = 1 if /All done/;
      push @buffer, $_;
      shift @buffer if scalar @buffer > 5;
    }
    close CMD;
  
    die "ERROR: Encountered an error running convert script\n".join("", @buffer)."\n" unless $finished;
  }
  
  $self->tar($type, $mod);
  
  return;
}


1;
