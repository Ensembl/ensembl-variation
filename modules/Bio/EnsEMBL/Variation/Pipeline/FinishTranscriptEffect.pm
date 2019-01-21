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
use FileHandle;

use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);

sub run {
  my $self = shift;

  my @transcript_ids = @{$self->param('transcript_ids_funnel')};

  my $var_dba = $self->get_species_adaptor('variation');
  my $tva = $var_dba->get_TranscriptVariationAdaptor;

 # initialise a hash of files
  my $files = {
    transcript_variation      => { 'cols' => [$tva->_write_columns],      },
    MTMP_transcript_variation => { 'cols' => [$tva->_mtmp_write_columns], },
  };

  # create filenames, file handles etc
  my $tmpdir = $self->param('pipeline_dir');
  $ImportUtils::TMP_DIR = $tmpdir;
  # COPY FILES FIRST IN CASE OF ERROR?
  # create file handles
  foreach my $transcript_id (@transcript_ids) {
    for my $table(keys %$files) {
      my $hash = $files->{$table};
      $hash->{filename} = sprintf('%s_%s.txt', $transcript_id, $table);
      $hash->{filepath} = sprintf('%s/%s', $tmpdir, $hash->{filename});
    }
    foreach my $table(keys %$files) {
      $ImportUtils::TMP_FILE = $files->{$table}->{filename};
      load($var_dba->dbc, ($table, @{$files->{$table}->{cols}}));
    }
  }

  return;
}

1;

