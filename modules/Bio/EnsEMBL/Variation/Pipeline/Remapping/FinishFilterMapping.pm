#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.




=head1 CONTACT

Please email comments or questions to the public Ensembl
developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

Questions may also be sent to the Ensembl help desk at
<helpdesk.org>.

=cut
package Bio::EnsEMBL::Variation::Pipeline::Remapping::FinishFilterMapping;

use strict;
use warnings;
use FileHandle;

use base ('Bio::EnsEMBL::Hive::Process');

sub fetch_input {
  my $self = shift;
}

sub run {
  my $self = shift;

  my $working_dir = $self->param('pipeline_dir');
  my $statistics_dir = "$working_dir/statistics";
  my $overall_counts = {};
  opendir(DIR, $statistics_dir) or die $!;
  while (my $file = readdir(DIR)) {
    if ($file =~ m/\.txt$/) {
      my $fh = FileHandle->new("$statistics_dir/$file", 'r');
      while (<$fh>) {
        chomp;
        my ($stats, $count) = split/=/;
        $overall_counts->{$stats} += $count;
      }
      $fh->close();
    }
  }
  closedir(DIR); 

  my $fh = FileHandle->new("$working_dir/overall_counts.txt", 'w');
  while (my ($stats, $counts) = each %$overall_counts) {
    print $fh "$stats=$counts\n";
  }
  $fh->close();
}


1;
