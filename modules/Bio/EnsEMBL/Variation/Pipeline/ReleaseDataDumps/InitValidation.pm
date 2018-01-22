=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

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
package Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::InitValidation;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::BaseDataDumpsProcess');

sub run {
  my $self = shift;
  my $file_type    = $self->param('file_type');
  my $species      = $self->param('species');
  my $pipeline_dir = $self->data_dir($species);

  my $working_dir = "$pipeline_dir/$file_type/$species/";
  my $files = get_files($working_dir);
  my @input = ();
  foreach my $file_name (keys %$files) {
    my $params = {};
    $params->{file_name}   = $file_name;
    $params->{working_dir} = $working_dir;
    push @input, $params;
  }
  $self->param('input_for_validation', \@input);
}

sub get_files {
  my $working_dir = shift;
  my $files = {};
  opendir(my $dh, $working_dir) or die $!;
  my @dir_content = readdir($dh);
  closedir($dh);
  foreach my $file (@dir_content) {
    if ($file =~ m/\.gvf$/) {
      $file =~ s/\.gvf//g;
      $files->{$file} = 1;    
    }
  }
  return $files;
}

sub write_output {
  my $self = shift;
  $self->dataflow_output_id($self->param('input_for_validation'), 2);
  return;
}


1;
