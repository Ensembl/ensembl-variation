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
package Bio::EnsEMBL::Variation::Pipeline::Remapping::BaseRemapping;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Hive::Process);

sub read_line {
  my $self = shift;
  my $line = shift;
  my @key_values = split("\t", $line);
  my $mapping = {};
  foreach my $key_value (@key_values) {
    my ($table_name, $value) = split('=', $key_value, 2);
    $mapping->{$table_name} = $value;
  }
  return $mapping;
}

sub count_files {
    my $self = shift;
    my $dir = shift;
    my $file_type = shift;
    opendir(my $dh, $dir) or die "Not a directory $dir";
    my $count = scalar(grep { $_ =~ m/$file_type$/ } readdir($dh));
    closedir($dh);
    return $count;
}

sub run_cmd {
    my $self = shift;
    my $cmd = shift;
    if (my $return_value = system($cmd)) {
        $return_value >>= 8;
        die "system($cmd) failed: $return_value";
    }
}

1;
