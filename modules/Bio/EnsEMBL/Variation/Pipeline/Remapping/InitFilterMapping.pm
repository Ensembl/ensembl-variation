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

use File::Path qw(make_path);

sub write_output {
    my $self = shift;

    my $mapping_results_dir = $self->param('mapping_results_dir');
    my $filtered_mappings_dir = $self->param('filtered_mappings_dir');
    my $load_features_dir = $self->param('load_features_dir');  
    if ($self->param('mode') eq 'remap_read_coverage') {
        my @input = ();
        opendir (IND_DIR, $mapping_results_dir) or die $!;
        while (my $individual_dir = readdir(IND_DIR)) {
            next if ($individual_dir =~ /^\./);
            make_path("$filtered_mappings_dir/$individual_dir") or die "Failed to create $filtered_mappings_dir/$individual_dir $!";
            make_path("$load_features_dir/$individual_dir") or die "Failed to create $load_features_dir/$individual_dir $!";

            my $count_mappings_files =  $self->compare_file_counts("$mapping_results_dir/$individual_dir");
            my $file_number = 1;
            while ($file_number <= $count_mappings_files) {
                die "No such file $mapping_results_dir/$individual_dir/mappings_$file_number.txt" unless (-e "$mapping_results_dir/$individual_dir/mappings_$file_number.txt");
                die "No such file $mapping_results_dir/$individual_dir/failed_mapping_$file_number.txt" unless (-e "$mapping_results_dir/$individual_dir/failed_mapping_$file_number.txt");
                push @input, {
                    'file_number' => $file_number,
                    'individual_id' => $individual_dir,
                };
                $file_number++;
            }
        }
        closedir(IND_DIR);
        $self->dataflow_output_id(\@input, 2);
    } else {
        my $count_mappings_files = $self->compare_file_counts($mapping_results_dir);
        my $input = $self->input_for_filter_mappings($mapping_results_dir, $count_mappings_files);
        $self->dataflow_output_id($input, 2);
    }
    1;
}

sub compare_file_counts {
    my $self = shift;
    my $dir = shift;
    my $count_mappings_files = 0;
    my $count_failed_mappings_files = 0;

    opendir(DIR, $dir) or die "Not a directory $dir";
    while (my $file = readdir(DIR)) {
        if ($file =~ m/^mappings_(.+)\.txt$/) {
            $count_mappings_files++;
        } elsif ($file =~ m/^failed_mapping_(.+)\.txt$/) {
            $count_failed_mappings_files++;
        }
    }
    closedir(DIR); 

    die "File count for mappings ($count_mappings_files) and failed_mappings ($count_failed_mappings_files) in $dir differs" unless ($count_mappings_files == $count_failed_mappings_files);
    return $count_mappings_files;
}

sub input_for_filter_mappigs {
    my $self = shift;
    my $dir = shift;
    my $count_mappings_files = shift;
    my @input;
    my $file_number = 1;

    while ($file_number <= $count_mappings_files) {
        die "No such file $dir/mappings_$file_number.txt" unless (-e "$dir/mappings_$file_number.txt");
        die "No such file $dir/failed_mapping_$file_number.txt" unless (-e "$dir/failed_mapping_$file_number.txt");
        push @input, {
            'file_number' => $file_number,
        };
        $file_number++;
    }
    return \@input; 
}



1;
