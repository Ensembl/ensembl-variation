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
package Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::CleanUp;

use strict;

use base ('Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::BaseDataDumpsProcess');

sub run {
  my $self = shift;
  my $tmp_dir = $self->param('tmp_dir');
  my $file_type = $self->param('file_type');

  my $species = $self->param('species');
  my $mode = $self->param('mode');
  my $data_dump_dir = $self->data_dir($species);

  if ($mode eq 'post_gvf_dump') {
    my $working_dir = "$data_dump_dir/$file_type/$species";
    opendir(my $dh, $working_dir) or die $!;
    my @dir_content = readdir($dh);
    closedir($dh);
    foreach my $file (@dir_content) {
      if ($file =~ m/gvf$/) {
        `gzip $working_dir/$file`;
      } 
      if ($file =~ m/^Validate/) {
        `mv $working_dir/$file $tmp_dir`;
      }
    }				
  }

  if ($mode eq 'post_join_dumps') {
    foreach my $file_type (qw/vcf gvf/) {
      my $working_dir = "$data_dump_dir/$file_type/$species";
      opendir(my $dh, $working_dir) or die $!;
      my @dir_content = readdir($dh);
      closedir($dh);
      foreach my $file (@dir_content) {
        if ($file =~ m/generic/) {
          my $file_name = $file;
          $file_name =~ s/_generic//;
          `mv $working_dir/$file $working_dir/$file_name`;
          if ($file_name =~ /gvf$/) {
            `gzip $working_dir/$file_name`;
          }
        } 
        if ($file =~ /gvf$/) {
          `gzip $working_dir/$file`;
        }
      }				
    }
  }
}

1;

