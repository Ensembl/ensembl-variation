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
package Bio::EnsEMBL::Variation::Pipeline::SpliceAI::BaseSpliceAI;

use strict;

use base ('Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess');

use Bio::EnsEMBL::Registry;
use File::Path qw(make_path);

sub checksum {
  my $path = shift;
  my $checksum = `sum $path`;
  $checksum =~ s/\s* $path//xms;
  chomp($checksum);
  return $checksum;
}

sub create_dir {
  my ($self, $dir) = @_;
  if (-d "$dir") {
    unless ($self->is_empty($dir)) {
      die("$dir is not empty. Delete files before running the pipeline.");
    }
  } else {
    make_path("$dir") or die "Failed to create dir $dir $!";
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

sub set_chr_from_filename {
  my $self = shift;
  my $vcf_file = $self->param('input_file');
  $vcf_file =~ /.*chr(.*).vcf/;
  my $chr = $1;
  $chr =~ s/\..*//;
  if (!$chr) {
    die("Could not get chromosome name from file name ($vcf_file).");
  }
  $self->param('chr', $chr);
}

1;

