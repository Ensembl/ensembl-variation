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
package Bio::EnsEMBL::Variation::Pipeline::DumpVEP::DumpVEP;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Pipeline::DumpVEP::BaseVEP);

sub param_defaults {
  return {
    'species_refseq' => 0,
    'species_flags'  => {},
    'overwrite'      => 0,
  };
}

sub run {
  my $self = shift;
  
  # basic params
  my $species  = $self->required_param('species');
  my $assembly = $self->required_param('assembly');
  my $version  = $self->required_param('ensembl_release');
  my $refseq   = $self->required_param('species_refseq') ? '--refseq' : '';
  my $dir      = $self->required_param('pipeline_dir');
  
  # species-specific
  my $species_flags = $self->param('species_flags');
  
  my $species_flags_cmd = $refseq.' ';
  if(my $flags = $species_flags->{$species}) {
    $species_flags_cmd .= join(' ', map {$flags->{$_} eq '1' ? '--'.$_ : '--'.$_.' '.$flags->{$_}} keys %$flags);
  }
  
  # db params
  my $host = $self->required_param('host');
  my $port = $self->required_param('port');
  my $user = $self->required_param('user');
  my $pass = $self->required_param('pass') ? '--pass '.$self->required_param('pass') : '';
  
  # cmnd line
  my $perl    = $self->required_param('perl_command');
  my $vep_dir = $self->required_param('ensembl_cvs_root_dir').'/ensembl-tools/scripts/variant_effect_predictor';
  my $vep     = $self->required_param('vep_command');
  
  # construct command
  my $cmd = sprintf(
    '%s %s/variant_effect_predictor.pl %s --host %s --port %i --user %s %s --species %s --assembly %s --db_version %s --dir %s %s',
    $perl,
    $vep_dir,
    $vep,
    
    $host,
    $port,
    $user,
    $pass,
    
    $species,
    $assembly,
    $version,
    $dir,
    $species_flags_cmd
  );
  
  my $finished = 0;
  
  open CMD, "$cmd 2>&1 |" or die "ERROR: Failed to run command $cmd";
  my @buffer;
  while(<CMD>) {
    $finished = 1 if /Finished/;
    push @buffer, $_;
    shift @buffer if scalar @buffer > 5;
  }
  close CMD;
  
  die "ERROR: Encountered an error running VEP\n".join("", @buffer)."\n" unless $finished;
  
  $self->tar($self->param('species_refseq') ? 'refseq' : '');
  
  return;
}


1;
