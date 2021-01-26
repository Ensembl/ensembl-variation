=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2021] EMBL-European Bioinformatics Institute

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
package Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::BaseDataDumpsProcess;

use strict;

use base ('Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess');

use Bio::EnsEMBL::Registry;
use File::Path qw(make_path);

sub compute_checksums_for_directory {
  my ($self, $working_dir) = @_;
  opendir(my $dh, $working_dir) or die "Cannot open directory $working_dir";
  my @files = sort {$a cmp $b} readdir($dh);
  closedir($dh) or die "Cannot close directory $working_dir";
  my @checksums = ();
  foreach my $file (@files) {
    next if $file =~ /^\./;
    next if $file =~ /^CHECKSUM/;
    my $path = File::Spec->catfile($working_dir, $file);
    my $checksum = checksum($path);
    push(@checksums, [$checksum, $file]);
  }
  my $fh = FileHandle->new("$working_dir/CHECKSUMS", 'w');
  foreach my $entry (@checksums) {
    my $line = join("\t", @{$entry});
    print $fh $line, "\n";
  }
  $fh->close();
}

sub checksum {
  my $path = shift;
  my $checksum = `sum $path`;
  $checksum =~ s/\s* $path//xms;
  chomp($checksum);
  return $checksum;
}

sub data_dir {
  my ($self,$species) = @_;
  my $data_dump_dir = $self->param('pipeline_dir');
  my $species_division = $self->param('species_division');
  # If division is defined append the pipeline_dir
  if ($species_division)
  {
    $data_dump_dir = $data_dump_dir."/".$species_division."/variation";
  }
  return $data_dump_dir;
}

sub create_species_dir {
  my ($self, $species_dir) = @_;
  if (-d "$species_dir") {
    unless ($self->is_empty($species_dir)) {
      die("$species_dir is not empty. Delete files before running the pipeline.");
    }
  } else {
    make_path("$species_dir") or die "Failed to create dir $species_dir $!";
  }
}

sub is_empty {
  my ($self, $dir) = @_;
  opendir(my $dh, $dir) or die "Not a directory: $dir.";
  my $count =  scalar(grep { $_ ne "." && $_ ne ".." } readdir($dh)) == 0;
  closedir($dh);
  return $count;
}

sub run_system_command {
  my ($self, $cmd) = @_; 
  my ($return_value, $stderr, $flat_cmd) = $self->SUPER::run_system_command($cmd);
  if ($return_value) {
    die("there was an error running as ($flat_cmd: $stderr)");
  }
}


1;

