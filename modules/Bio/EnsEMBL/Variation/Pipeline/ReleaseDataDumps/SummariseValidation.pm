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
package Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::SummariseValidation;

use strict;

use base ('Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::BaseDataDumpsProcess');

sub fetch_input {
  my $self = shift;
  my @input;
  my $file_type   = $self->param('file_type');
  my $species = $self->param('species');
  my $pipeline_dir = $self->data_dir($species);

  my $summary_file_name = '';
  if ($file_type eq 'gvf') {
    $summary_file_name = 'SummaryValidateGVF.txt';
  } elsif ($file_type eq 'vcf') {
    $summary_file_name = 'SummaryValidateVCF.txt';
  } else {
    die "File type ($file_type) is not correct. It must be gvf of vcf";
  }

  open(my $fh, '>>', "$pipeline_dir/$summary_file_name") or die "Could not open file '$pipeline_dir/$summary_file_name' $!";

  my $failed_dumps = 0;

  my $path = "$pipeline_dir/$file_type/$species";
  opendir(my $dh, $path) or die $!;
  my @dir_content = readdir($dh);
  closedir($dh);
  foreach my $file (@dir_content) {
    if ($file =~ m/^Validate(.)*out$/) {
      if ($file_type eq 'gvf') {
        my $return_value = `grep 'No Errors found in this file' $path/$file`;
        if (length($return_value) == 0) {
          print $fh "GVF validator reports errors for $file\n";
          $failed_dumps++;
        }
      } elsif ($file_type eq 'vcf') {
      }
    }
  }
  close $fh;
  $self->warning("For $species, $failed_dumps dumps report errors.") if ($failed_dumps > 0);
}

1;
