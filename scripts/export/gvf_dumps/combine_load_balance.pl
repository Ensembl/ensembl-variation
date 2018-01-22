# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute
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

use warnings;
use strict;

use FileHandle;

my $gvf_dir = '/dir_to_gvfs/';

my @file_names = (
    'Homo_sapiens_incl_consequences',
    '1000GENOMES-phase_1_AFR',
    '1000GENOMES-phase_1_ASN',
    '1000GENOMES-phase_1_AMR',
    '1000GENOMES-phase_1_EUR',
);

foreach my $file_name (@file_names) {

    my $count = 1;

    my $fh   = FileHandle->new($gvf_dir . $file_name . '.gvf', 'w') ;
    my $fh_1 = FileHandle->new($gvf_dir . 'completed/' . $file_name . ".1.gvf", 'r');
    while (<$fh_1>) {
        chomp;
        if (/^##/) {
            print $fh $_, "\n";
        } else {
            last;
        }
    }

    foreach my $i (2..12) {
        my $fh_tmp = FileHandle->new($gvf_dir . 'completed/' . $file_name . ".$i.gvf", 'r');
        while (<$fh_tmp>) {
            chomp;
            if (/^##/) {
                if (/^##sequence-region/) {
                    print $fh $_, "\n";
                }
            }  else {
                last;
            }
        }
    }

    foreach my $i (1..12) {
        my $fh_tmp = FileHandle->new($gvf_dir . 'completed/' . $file_name . ".$i.gvf", 'r');
        while (<$fh_tmp>) {
            chomp;
            next if /^##/;
            my @feature_lines = split(/\s+/, $_, 9);
            my $attributes = $feature_lines[8];
            my @attributes_split = split(/;/, $attributes, 2);
            print $fh join("\t", @feature_lines[0..7]), "\tID=$count;", $attributes_split[1], "\n";
            $count++;
        }
    } 
}

__END__

=head1 NAME

combine_load_balance.pl

=head1 DESCRIPTION

Combine gvf snippets.

=head1 SYNOPSIS

combine_load_balance.pl --gvf_output FILE --header_file FILE --log_file FILE

=head1 EXAMPLE COMMAND LINES

=head1 OPTIONS

=over 4

=item B<--gvf_output FILE>

=item B<--header_file FILE>

=item B<--log_file FILE>

=item B<--help>

Display this documentation

=head1

For help with this script address questions to http://lists.ensembl.org/mailman/listinfo/dev
