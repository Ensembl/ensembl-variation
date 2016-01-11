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
package Bio::EnsEMBL::Variation::Pipeline::FinishTranscriptEffect;

use strict;
use warnings;
use ImportUtils qw(load);
use Sys::Hostname;

use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);

my $DEBUG = 0;

sub run {
  my $self = shift;

  my $dbc = $self->get_species_adaptor('variation')->dbc;

  my $dir = $self->required_param('pipeline_dir');
  $ImportUtils::TMP_DIR = $dir;

  my $host = hostname;

  # do unique sort on command line, it's faster than relying on MySQL's unique index
  foreach my $file(grep {-e "$dir/$_"} qw(variation_hgvs.txt variation_genename.txt)) {
    system(
      sprintf(
        'sort %s -u %s/%s > %s/%s.unique',
        $host =~ /sanger/ ? '--parallel=4' : '', $dir, $file, $dir, $file
      )
    ) and die("ERROR: Failed to unique sort $file");
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

1;

