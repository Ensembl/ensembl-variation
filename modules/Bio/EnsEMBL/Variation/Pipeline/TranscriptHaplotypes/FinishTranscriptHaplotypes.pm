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

package Bio::EnsEMBL::Variation::Pipeline::TranscriptHaplotypes::FinishTranscriptHaplotypes;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);

use FileHandle;

my $DEBUG = 0;

sub run {
  my $self = shift;

  my $dir = $self->required_param('pipeline_dir');

  my $handles = {};

  for my $t1(qw(cds protein)) {
    for my $t2(qw(hap dip)) {
      my $fh = FileHandle->new();
      $fh->open(">$dir/".$t1."_".$t2."lotypes.txt") or die $!;
      $handles->{$t1.'_'.$t2} = $fh;
    }
  }

  opendir DIR, $dir.'/haplotypes/';

  foreach my $hex_stub(grep {!/^\./} readdir DIR) {

    opendir HEX, "$dir/haplotypes/$hex_stub";

    my @files = grep {!/^\./} readdir HEX;

    foreach my $handle_key(keys %$handles) {
      my $fh = $handles->{$handle_key};
      foreach my $file(grep {$_ =~ /$handle_key/} @files) {
        open IN, "$dir/haplotypes/$hex_stub/$file" or die $!;
        while(<IN>) {
          print $fh $_;
        }
        close IN;
      }
    }
  }
}

1;
