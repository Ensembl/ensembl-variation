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
package Bio::EnsEMBL::Variation::Pipeline::DumpVEP::BaseVEP;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);

sub tar {
  my $self = shift;
  my $type = shift;
  my $mod  = shift;
  
  my $debug    = $self->param('debug');
  my $species  = $self->required_param('species');
  my $assembly = $self->required_param('assembly');
  my $version  = $self->required_param('ensembl_release');
  my $dir      = $self->required_param('pipeline_dir');
  
  $species .= $type ? '_'.$type : '';
  $mod ||= '';
  
  my $tar_file = sprintf(
    '%s/%s_vep_%i_%s%s.tar.gz',
    $dir,
    $species,
    $version,
    $assembly,
    $mod
  );
  
  # check if tar exists
  if(!$self->param('overwrite') && -e $tar_file) {
    print STDERR "Existing dump file found for $species, skipping (use --overwrite to overwrite)\n";
    return;
  }
  
  # check dir exists
  my $root_dir = $dir;
  my $sub_dir  = $species."/".$version."_".$assembly;
  
  die("ERROR: VEP dump directory $root_dir/$sub_dir not found") unless -e $root_dir.'/'.$sub_dir;
  
  my $command = "tar -cz -C $root_dir -f $tar_file $sub_dir";
  
  if($debug) {
    print STDERR "$command\n";
  }
  else {
    my $output = `$command`;
    die "ERROR: Failed to create tar file $tar_file\n$output\n" if $output;
  }
  
  return;
}

1;
