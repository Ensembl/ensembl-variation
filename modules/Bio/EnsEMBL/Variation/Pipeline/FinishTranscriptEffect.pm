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
package Bio::EnsEMBL::Variation::Pipeline::FinishTranscriptEffect;

use strict;
use warnings;
use ImportUtils qw(load);
use FileHandle;

use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);

sub run {
  my $self = shift;

  my @transcript_ids = @{$self->param('transcript_ids_funnel')};
  my $gene_stable_id = $self->param('gene_stable_id');

  my $var_dba = $self->get_species_adaptor('variation');
  my $tva = $var_dba->get_TranscriptVariationAdaptor;

  my $files = {
    transcript_variation      => { 'cols' => [$tva->_write_columns],      },
    MTMP_transcript_variation => { 'cols' => [$tva->_mtmp_write_columns], },
  };

  # create log file and read if available or create

  my $pipeline_dir = $self->param('pipeline_dir');
  my $loaded_transcript_ids = {};
  my $fh;
  my $log_file = "$pipeline_dir/load_log_files/LOG_$gene_stable_id";
  if (-e $log_file) {
    # read successfully loaded transcripts
    $fh = FileHandle->new($log_file, 'r');
    while (<$fh>) {
      chomp;
      $loaded_transcript_ids->{$_} = 1;
    }
    $fh->close;
  }
  open($fh, '>>', $log_file) or die "Could not open file '$log_file' $!";

  foreach my $transcript_id (@transcript_ids) {
    foreach my $table(keys %$files) {
      next if $loaded_transcript_ids->{"$transcript_id\_$table"};
      my $tmpdir =  $self->get_files_dir($transcript_id, 'transcript_effect');
      my $filename = sprintf('%s_%s.txt', $transcript_id, $table);
      $self->run_cmd("cp $tmpdir/$filename $tmpdir/load_$filename");      
      $ImportUtils::TMP_DIR = $tmpdir;
      $ImportUtils::TMP_FILE = "load_$filename";
      load($var_dba->dbc, ($table, @{$files->{$table}->{cols}}));
      print $fh "$transcript_id\_$table\n";
      $self->run_cmd("rm $tmpdir/$filename");
    }
  }

  close($fh);

  return;
}

1;

