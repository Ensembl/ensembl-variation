=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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
package Bio::EnsEMBL::Variation::Pipeline::DumpVEP::DistributeDumps;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Pipeline::DumpVEP::BaseVEP);

use FileHandle;
use File::Path qw(make_path);
use File::Spec;

sub run {
  my $self = shift;
  my $version = $self->required_param('ensembl_release');
  my $dir = $self->required_param('pipeline_dir');

  foreach my $folder (qw/web production rest/) {
    make_path("$dir/$folder") if (!-d "$dir/$folder");
  }

  # rest
  foreach my $assembly (qw/GRCh37 GRCh38/) {
    my $file = "homo_sapiens_vep_$version\_$assembly\_tabixconverted.tar.gz";
    if (-e "$dir/$file") {
      run_cmd("ln $dir/$file $dir/rest");
    } else {
      $self->warning("$file doesn't exist");
    }
  }

  # production web 
  opendir (DIR, $dir) or die $!;
  while (my $file = readdir(DIR)) {
    if ($file =~ /gz$/ && $file !~ /tabix/) {
      run_cmd("ln $dir/$file $dir/production/");
      run_cmd("ln $dir/$file $dir/web/");
    }
  }
  closedir (DIR) or die $!;

  opendir (DIR, $dir) or die $!;
  while (my $file = readdir(DIR)) {
    if ($file =~ /tabix/) {
      my $file_no_tabix = $file;
      $file_no_tabix =~ s/\_tabixconverted//;
      run_cmd("rm $dir/web/$file_no_tabix");  
      run_cmd("ln $dir/$file $dir/web/$file_no_tabix");  
    }
  }
  closedir (DIR) or die $!;

  compute_checksums("$dir/production/");

}

sub compute_checksums {
  my $dir = shift;
  opendir(my $dh, $dir) or die $!;
  my @files = sort {$a cmp $b} readdir($dh);
  closedir($dh) or die $!;
  my @checksums = ();
  foreach my $file (@files) {
    next if $file =~ /^\./;
    next if $file =~ /^CHECKSUM/;
    my $path = File::Spec->catfile($dir, $file);
    my $checksum = checksum($path);
    push(@checksums, [$checksum, $file]);
  }
  my $fh = FileHandle->new("$dir/CHECKSUMS", 'w');
  foreach my $entry (@checksums) {
    my $line = join("\t", @{$entry});
    print $fh $line, "\n";
  }
  $fh->close();
}

sub run_cmd {
  my $cmd = shift;
  if (my $return_value = system($cmd)) {
    $return_value >>= 8;
    die "system($cmd) failed: $return_value";
  }
}

sub checksum {
  my $path = shift;
  my $checksum = `sum $path`;
  $checksum =~ s/\s* $path//xms;
  chomp($checksum);
  return $checksum;
}

1;
