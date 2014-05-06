#!/usr/bin/env perl
# Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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
package Bio::EnsEMBL::Variation::Pipeline::Remapping::InitFilterMapping;

use strict;

use base ('Bio::EnsEMBL::Hive::Process');

sub write_output {
    my $self = shift;

    my $mapping_results_dir = $self->param('mapping_results_dir');
    
    my $count_mappings_files = 0;
    my $count_failed_mappings_files = 0;

    opendir(DIR, $mapping_results_dir) or die "Not a directory $mapping_results_dir";
    while (my $file = readdir(DIR)) {
        if ($file =~ m/^mappings_(.+)\.txt$/) {
            $count_mappings_files++;
        } elsif ($file =~ m/^failed_mapping_(.+)\.txt$/) {
            $count_failed_mappings_files++;
        }
    }
    closedir(DIR); 

    die "File count for mappings ($count_mappings_files) and failed_mappings ($count_failed_mappings_files) in $mapping_results_dir differs" unless ($count_mappings_files == $count_failed_mappings_files);

    my @input;
    my $file_number = 1;

    while ($file_number <= $count_mappings_files) {
        die "No such file $mapping_results_dir/mappings_$file_number.txt" unless (-e "$mapping_results_dir/mappings_$file_number.txt");
        die "No such file $mapping_results_dir/failed_mapping_$file_number.txt" unless (-e "$mapping_results_dir/failed_mapping_$file_number.txt");
        push @input, {
            'file_number' => $file_number,
        };
        $file_number++;
    }
    $self->dataflow_output_id(\@input, 2);
    1;
}


1;
