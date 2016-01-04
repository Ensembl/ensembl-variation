=head1 LICENSE

Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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
package Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::Finish;

use strict;
use warnings;
use base ('Bio::EnsEMBL::Variation::Pipeline::ReleaseDataDumps::BaseDataDumpsProcess');

use FileHandle;
use File::Path qw(make_path);
use File::Spec;

sub fetch_input {
  my $self = shift;
}

sub run {
  my $self = shift;
  my $data_dump_dir = $self->param('pipeline_dir');
  my $all_species = $self->get_all_species();
  $self->delete_tmp_files($data_dump_dir, $all_species);
  tabix_vcf_files($data_dump_dir, $all_species);
  $self->add_readme($data_dump_dir, $all_species);
  compute_checksums($data_dump_dir, $all_species);
  lc_dir($data_dump_dir, $all_species);
  rm_uc_dir($data_dump_dir, $all_species);
}

sub delete_tmp_files {
  my ($self, $data_dir, $all_species) = @_;
  my $tmp_dir = $self->param('tmp_dir');
  foreach my $species (keys %$all_species) {
    foreach my $file_type (qw/vcf gvf/) {
      my $working_dir = "$data_dir/$file_type/$species/"; 
      opendir(DIR, $working_dir) or die $!;
      while (my $file = readdir(DIR)) {
        if ($file =~ m/_validate\.vcf.gz$|\.txt$/) {
          system("mv $working_dir/$file $tmp_dir");
          $self->warning($file); 
        }
      }
      closedir(DIR);    
    }
  }
}

sub tabix_vcf_files {
  my ($data_dir, $all_species) = @_;
  foreach my $species (keys %$all_species) {
    my $working_dir = "$data_dir/vcf/$species/"; 
    opendir(DIR, $working_dir) or die $!;
    while (my $file = readdir(DIR)) {
      if ($file =~ m/\.vcf.gz$/) {
        system("tabix -p vcf $working_dir/$file");
      }
    }
    closedir(DIR);    
  }
}

sub add_readme {
  my ($self, $data_dir, $all_species) = @_;
  foreach my $file_type (qw/gvf vcf/) {
    my $readme_file = $self->param("$file_type\_readme");
    my $fh = FileHandle->new($readme_file, 'r');
    my $readme;
    {
      local $/ = undef;
      $readme = <$fh>;
    }    
    $fh->close();    
    foreach my $species (keys %$all_species) {
      open README, ">$data_dir/$file_type/$species/README" or die "Failed to create README file for species $species\n";
      print README sprintf($readme, $species, $species, $species, $species), "\n";
      close README;
    }
  }
}

sub compute_checksums {
  my ($data_dir, $all_species) = @_;  

  foreach my $file_type (qw/vcf gvf/) {
    foreach my $species (keys %$all_species) {
      my $working_dir = "$data_dir/$file_type/$species/";
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
  }
}

sub checksum {
  my $path = shift;
  my $checksum = `sum $path`;
  $checksum =~ s/\s* $path//xms;
  chomp($checksum);
  return $checksum;
}

sub lc_dir {
  my ($data_dir, $all_species) = @_;
  foreach my $species (keys %$all_species) {
    my $lc_species = lc $species;
    foreach my $file_type (qw/gvf vcf/) {
      my $lc_dir = "$data_dir/$file_type/$lc_species/"; 
      make_path($lc_dir) unless (-d $lc_dir);
      system("mv $data_dir/$file_type/$species/* $lc_dir");
    }
  }
}

sub rm_uc_dir {
  my ($data_dir, $all_species) = @_;
  foreach my $species (keys %$all_species) {
    foreach my $file_type (qw/gvf vcf/) {
      my $dir = "$data_dir/$file_type/$species/";
      system("rm -r $dir");
    }
  }
}

1;
