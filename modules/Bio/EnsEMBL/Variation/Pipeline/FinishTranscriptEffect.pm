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
package Bio::EnsEMBL::Variation::Pipeline::FinishTranscriptEffect;

use strict;
use warnings;
use ImportUtils qw(load);
use Sys::Hostname;
use FileHandle;
use File::Path qw(rmtree);

use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);

my $DEBUG = 0;

sub run {
  my $self = shift;

  $self->rejoin_table_files();

  my $dbc = $self->get_species_adaptor('variation')->dbc;

  my $dir = $self->required_param('pipeline_dir');
  $ImportUtils::TMP_DIR = $dir;

  my $host = hostname;

  # do unique sort on command line, it's faster than relying on MySQL's unique index
  foreach my $file(grep {-e "$dir/$_"} qw(variation_hgvs.txt variation_genename.txt)) {
    system("gzip -c $dir/$file > $dir/$file\_bak.gz");
    system(
      sprintf(
        'cat %s/%s | sort -T %s -u > %s/%s.unique',
        $dir, $file, $dir, $dir, $file, 
      )
    ) and die("ERROR: Failed to unique sort $file");
    unlink("$dir/$file\.gz") if -e "$dir/$file\.gz";
    system("gzip $dir/$file");# unlink("$dir/$file");
  }

  if(-e $dir.'/variation_hgvs.txt.unique') {
    $ImportUtils::TMP_FILE = 'variation_hgvs.txt.unique';
    load($dbc, qw(variation_hgvs variation_id hgvs_name));
  }

  if(-e $dir.'/variation_genename.txt.unique') {
    $ImportUtils::TMP_FILE = 'variation_genename.txt.unique';
    load($dbc, qw(variation_genename variation_id gene_name));
  }

  return;
}

sub rejoin_table_files {
  my $self = shift;

  my $dir = $self->required_param('pipeline_dir');

  my $gene_fh = FileHandle->new();
  $gene_fh->open(">".$dir."/variation_genename.txt") or die $!;
  my $hgvs_fh = FileHandle->new();
  $hgvs_fh->open(">".$dir."/variation_hgvs.txt") or die $!;

  opendir DIR, $dir."/table_files";
  foreach my $hex_stub(grep {!/^\./} readdir DIR) {

    opendir HEX, "$dir/table_files/$hex_stub";
    foreach my $file(grep {!/^\./} readdir HEX) {

      my $fh = $file =~ /hgvs/ ? $hgvs_fh : $gene_fh;

      open IN, "$dir/table_files/$hex_stub/$file" or die $!;
      while(<IN>) {
        print $fh $_;
      }
      close IN;
    }
  }

  rmtree($dir."/table_files");
}

1;

