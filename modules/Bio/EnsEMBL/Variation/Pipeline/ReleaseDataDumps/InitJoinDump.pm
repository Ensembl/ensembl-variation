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
package Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::InitJoinDump;

use strict;
use warnings;

use FileHandle;

use base ('Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::BaseDataDumpsProcess');

sub run {
  my $self = shift;
  my $species      = $self->param('species');
  my $pipeline_dir = $self->data_dir($species);


  my @input = ();
 
  foreach my $file_type (qw/gvf vcf/) {

    my $dir = "$pipeline_dir/$file_type/$species/";
    my $files = {};

    opendir(my $dh, $dir) or die $!;
    my @dir_content = readdir($dh);
    closedir($dh);
    foreach my $file (@dir_content) {
      next if ($file =~ m/^\./);
      if ($file =~ m/\.$file_type\.gz/) {
        my $file_name = $file;
        $file_name =~ s/\.$file_type\.gz//;
        my ($dump_type, $range) = split('-', $file_name);
        push @{$files->{$dump_type}}, $range if ($range);
      }
    }

    foreach my $dump_type (keys %$files) {
      my $file_name = $dump_type;
      if ($dump_type =~ /_generic/) {
        $file_name =~ s/_generic//;
      }
      if (scalar @{$files->{$dump_type}} > 0) {
        push @input, {
          dir => $dir,
          file_name => $file_name,
          dump_type => $dump_type,
          input_ids => $files->{$dump_type},
          file_type => $file_type,
          mode => 'final_join',
        };
      }
    }
  }
  $self->param('input', \@input);
}


sub write_output {
  my $self = shift;
  if (scalar @{$self->param('input')} > 0) {
    $self->dataflow_output_id($self->param('input'), 2);
  } else {
    $self->dataflow_output_id([{mode => 'no_join'}], 2);
  }
  return;
}


1;
